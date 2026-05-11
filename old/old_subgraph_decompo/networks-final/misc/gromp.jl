using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf


#general
includet("../../utils/import_utils.jl")

# utils colge
includet("utils/utils-subgraphdecompo.jl")
includet("utils/partition-graph.jl")
includet("utils/checkers.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")

includet("../../heuristics/mepso.jl")
includet("../../compact/compact_plus.jl")

function solve_gromp(instance)
    
    
    # Budget : 60 seconds
    time_submappings = 30
    time_cg_heuristic = 30

    v_network = instance.v_network
    s_network = instance.s_network

    println("Starting...")
    time_beginning = time()



    # ======= SETTING UP THE DECOMPOSITION ======= #
    nb_virtual_subgraph = floor(Int, nv(v_network.graph)/10)
    v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.2)

    vn_decompo = set_up_decompo(instance, v_node_partitionning)
    
    println("Decomposition set: ")
        println("For $v_network, there is $(length(vn_decompo.subgraphs)) subgraphs:")

    for subgraph in vn_decompo.subgraphs
        println("       $(subgraph.graph[][:name]) with $(nv(subgraph.graph)) nodes")
    end
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")

    




    # ====== PAVING THE NETWORK WITH HEURISTIC ======= #

    println("Paving time...")
    time_0 = time()

    
    sub_mappings = find_submappings(instance, vn_decompo, solver="mepso")
    println("Mappings gotten! In just $(time() - time_0)")


    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    print("Master problem set... ")
    for v_subgraph in vn_decompo.subgraphs
        for mapping in sub_mappings[v_subgraph]
            add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
        end
    end
    print("Submappings added...")

    
    # ======= GETTING A SOLUTION ======= #
    value_cg_heuristic, cg_heuristic_solution = basic_heuristic(instance, vn_decompo, master_problem, time_cg_heuristic)


    result = Dict()
    result["solving_time"] = time() - time_beginning
    result["mapping_cost"] = value_cg_heuristic

    return result
end





function find_submappings(instance, vn_decompo; solver="mepso")


    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    size_max_v_subgraph = maximum(nv(v_subgraph.graph) for v_subgraph in vn_decompo.subgraphs)
    nb_cluster = floor(Int, nv(s_network) / (size_max_v_subgraph*1.5))

    mappings = []
    
    mappings_per_subgraph = Dict()
    for v_subgraph in vn_decompo.subgraphs
        mappings_per_subgraph[v_subgraph] = []
    end

    partitionning = partition_kahip(s_network.graph, nb_cluster, 0.15)
    clusters = []
    for i_subgraph in 1:nb_cluster
        cluster = Vector{Int64}() 
        for s_node in vertices(s_network)
            if partitionning[s_node] == i_subgraph
                push!(cluster, s_node)
            end
        end
        print("Cluster $i_subgraph has $(length(cluster)) nodes ")
        push!(clusters, cluster)
    end
    

    while length(mappings) < 200    
        for (i_cluster, cluster) in enumerate(clusters) 
            induced_subg = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")
            s_subgraph = Subgraph(induced_subg, cluster)

            for v_subgraph in vn_decompo.subgraphs 

                sub_instance = Instance(v_subgraph.graph, s_subgraph.graph)

                if solver=="mepso"
                    sub_mapping, cost = solve_mepso(sub_instance; nb_particle=30, nb_iter=50, time_max=0.25, print_things=false)
                elseif solver=="ilp"
                    result_milp = solve_compact_ffplus(sub_instance; time_solver = 30, stay_silent=true, linear=false)
                    sub_mapping = result_milp["mapping"]
                else
                    println("I don't know your solver, using mepso")
                    sub_mapping, cost = solve_mepso(sub_instance; nb_particle=30, nb_iter=50, time_max=0.25, print_things=false)
                end

                if isnothing(sub_mapping)
                    continue
                end

                true_cost = 0
                node_placement = []
                for v_node in vertices(v_subgraph.graph)
                    real_s_node = s_subgraph.nodes_of_main_graph[sub_mapping.node_placement][v_node]
                    append!(node_placement, real_s_node)
                    true_cost += s_network[real_s_node][:cost]
                end
    
    
                edge_routing = Dict()
                for v_edge in edges(v_subgraph.graph)
                    used_edges = []
                    for s_edge in sub_mapping.edge_routing[v_edge].edges
                        real_s_edge = get_edge(s_network_dir, s_subgraph.nodes_of_main_graph[src(s_edge)], s_subgraph.nodes_of_main_graph[dst(s_edge)])
                        push!(used_edges, real_s_edge)
                        true_cost += s_network_dir[src(real_s_edge), dst(real_s_edge)][:cost]
                    end
                    edge_routing[v_edge] = order_path(s_network_dir, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
                end
                real_sub_mapping = Mapping(v_subgraph.graph, s_network_dir, node_placement, edge_routing)
    
                push!(mappings, real_sub_mapping)
                push!(mappings_per_subgraph[v_subgraph], real_sub_mapping)

            end

        end

    end

    println("We have $(length(mappings)) mappings!")
    return mappings_per_subgraph    
end




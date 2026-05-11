
using Revise
using Statistics
using Graphs, MetaGraphsNext

includet("utils/utils-subgraphdecompo.jl")
includet("../../heuristics/uepso.jl")
includet("../../heuristics/mepso.jl")
includet("utils/partition-graph.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")



function solve_quitefast(instance) 



    # hum
    time_cg_heuristic = 5


    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    println("Starting...")
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network


    # ======= SETTING UP THE DECOMPOSITION ======= #
    nb_virtual_subgraph = floor(Int, nv(v_network.graph)/10)
    #v_node_partitionning = partition_graph_metis(instance.v_network.graph, nb_virtual_subgraph)
    v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.2)

    vn_decompo = set_up_decompo(instance, v_node_partitionning)
    
    println("Decomposition set: ")
        println("For $v_network, there is $(length(vn_decompo.subgraphs)) subgraphs:")

    for subgraph in vn_decompo.subgraphs
        println("       $(subgraph.graph[][:name]) with $(nv(subgraph.graph)) nodes")
    end
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")
    vn_subgraphs = vn_decompo.subgraphs

    
    # Set up master problem !
    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    print("Master problem set... ")


    # Decompose the substrate network into the same number of subgraph
    sn_subgraphs = []
    s_node_partitionning = partition_graph(instance.s_network.graph, nb_virtual_subgraph,  max_umbalance=1.3)
    for (i_cluster, cluster) in enumerate(s_node_partitionning)
        sub_s_network = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")
        push!(sn_subgraphs, Subgraph(sub_s_network, cluster))
    end

    # Map the v_subgraph into s_subgraph
    
    mappings_result = Dict()
    for v_subgraph in vn_subgraphs
        
        mappings = []

        for s_subgraph in sn_subgraphs
            sub_instance = Instance(v_subgraph.graph, s_subgraph.graph)
                
            sub_mapping, cost = solve_mepso(sub_instance; nb_particle=25, nb_iter=50, time_max=0.25, print_things=false)
            # Need to correct the mapping here
            
            if isnothing(sub_mapping) # invalid submapping!
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
            add_column(master_problem, instance, v_subgraph, real_sub_mapping, true_cost)

        end

        mappings_result[v_subgraph] = mappings


    end        

    # ======= GETTING A SOLUTION ======= #
    value_cg_heuristic, heuristic_mapping = basic_heuristic(instance, vn_decompo, master_problem, time_cg_heuristic)

    result = Dict()
    result["mapping"] = heuristic_mapping
    result["time_solving"] = time()-time_beginning
    result["value_cg_heuristic"] = round(Int, value_cg_heuristic)

    return result
end










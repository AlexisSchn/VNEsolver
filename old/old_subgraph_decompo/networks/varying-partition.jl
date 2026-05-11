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

function solve_varying(instance)
    
    
    # Budget : 60 seconds
    time_init = 30
    time_cg_heuristic = 30


    println("Starting...")
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network


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

    
    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    print("Master problem set... ")




    # ====== PAVING THE NETWORK WITH HEURISTIC ======= #

    println("Paving time...")
    time_0 = time()

    
    mappings = init_varying_partition(instance, vn_decompo)
    println("Mappings gotten! In just $(time() - time_0)")
    for v_subgraph in vn_decompo.subgraphs
        for mapping in mappings[v_subgraph]
            add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
        end
    end


        

    
    # ======= GETTING A SOLUTION ======= #
    value_cg_heuristic, cg_heuristic_solution = basic_heuristic(instance, vn_decompo, master_problem, time_cg_heuristic)
    #local_search(instance, vn_decompo, heur_sol)



    result = Dict()
    result["algo"] = "SHEEESH"
    result["time_solving"] = time() - time_beginning
    result["value_cg_heuristic"] = value_cg_heuristic

    return result
end





function init_varying_partition(instance, vn_decompo)


    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    size_max_v_subgraph = maximum(nv(v_subgraph.graph) for v_subgraph in vn_decompo.subgraphs)
    nb_cluster = floor(Int, nv(s_network) / (size_max_v_subgraph*1.5))

    mappings = []
    
    mappings_per_subgraph = Dict()
    for v_subgraph in vn_decompo.subgraphs
        mappings_per_subgraph[v_subgraph] = []
    end

    while length(mappings) < 300
        clusters = partition_graph_kahip_random_seed(s_network.graph, nb_cluster)
        for (i_cluster, cluster) in enumerate(clusters) 
            println("Cluster $i_cluster has $(length(cluster)) nodes ")
            induced_subg = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")
            s_subgraph = Subgraph(induced_subg, cluster)

            for v_subgraph in vn_decompo.subgraphs 

                sub_instance = Instance(v_subgraph.graph, s_subgraph.graph)
            

                sub_mapping, cost = solve_mepso(sub_instance; nb_particle=30, nb_iter=50, time_max=0.25, print_things=false)

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




function partition_graph_kahip_random_seed(graph, nb_clusters)

    
    # 1 : Partitionner
    inbalance = 0.1
    partition = disturbed_kahip_seed(graph, nb_clusters, inbalance)
    println("Partition of SN graph: $partition")
    clusters = [Vector{Int64}() for i in 1:nb_clusters]
    for s_node in vertices(graph)
        push!(clusters[partition[s_node]], s_node)
    end

    # 2 : Corriger
    for cluster in clusters
        simple_subgraph, vmap = induced_subgraph(graph, cluster)
        if !is_connected(simple_subgraph)
            #print("Issue with unconnected subgraphs to correct...")
            components = connected_components(simple_subgraph)
            component_sorted = sort(components, by=x->length(x), rev=true)
            for subcluster in component_sorted[2:length(component_sorted)]
                nodes_original = [vmap[node] for node in subcluster]
                subgraph_neighbors = zeros(Int, nb_clusters)
                for node in nodes_original
                    for neighbor in neighbors(graph, node)
                        if neighbor ∉ cluster
                            subgraph_neighbors[partition[neighbor]] += 1
                        end
                    end
                end
                most_connected_subgraph = sortperm(subgraph_neighbors, rev=true)
                cluster_to_put_nodes_in = most_connected_subgraph[1]
                append!(clusters[cluster_to_put_nodes_in], nodes_original)
                for node in nodes_original
                    partition[node] = cluster_to_put_nodes_in
                end
                filter!(e->e∉nodes_original, cluster)
            end
        end
    end
    
    

    return clusters

end


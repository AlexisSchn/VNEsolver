-656
using Revise
using Statistics
using Graphs, MetaGraphsNext

includet("utils/utils-subgraphdecompo.jl")

includet("utils/partition-graph.jl")

includet("heuristics/mepso-enhanced.jl")



function solve_sbire(instance) 

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    println("Starting...")
    time_beginning = time()


    # ======= SETTING UP THE DECOMPOSITION ======= #
    nb_virtual_subgraph = floor(Int, nv(v_network.graph)/10)
    v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.2)

    vn_decompo = set_up_decompo(instance, v_node_partitionning)
    
    print_stuff_decompo(vn_decompo, instance)
    
    vn_subgraphs = vn_decompo.subgraphs


    # shortest path things
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network_dir)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir, s_edge, :cost)
    end
    base_shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)
    

    # Decompose the substrate network into the same number of subgraph
    sn_subgraphs = []
    clusters = partition_graph(instance.s_network.graph, nb_virtual_subgraph, max_umbalance=1.25)
    for (i_cluster, cluster) in enumerate(clusters)
        sub_s_network = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")
        push!(sn_subgraphs, Subgraph(sub_s_network, cluster))
    end


    # Assign virtual subgraphs to substrate subgraphs
    assignment_virtual_substrate_subgraphs = Dict()
    # Large to large, in terms of density and number of nodes ?
    densities_v_subgraphs = [ne(v_subgraph.graph) for v_subgraph in vn_subgraphs]
    densities_s_subgraphs = [ne(s_subgraph.graph) for s_subgraph in sn_subgraphs]
    sorted_v_subgraphs = sortperm(densities_v_subgraphs)
    sorted_s_subgraphs = sortperm(densities_s_subgraphs)

    for i_rank in 1:(nb_virtual_subgraph) 
        # Weirdly, it doesnt help to consider the densities or other things. What works the best is the basic order...
        # SO WEIRD.
        v_subgraph = vn_subgraphs[i_rank]
        s_subgraph = sn_subgraphs[i_rank]
        assignment_virtual_substrate_subgraphs[v_subgraph] = s_subgraph
    end


    #print(assignment_virtual_substrate_subgraphs)
    
    # Get most central node of each substrate subgraph. 
    # Temporary placement: all virtual nodes are placed on the most central node of the substrate subgraph they are placed on.
    temporary_placement = zeros(Int, nv(v_network))
    for v_subgraph in vn_subgraphs
        s_subgraph = assignment_virtual_substrate_subgraphs[v_subgraph]
        cluster = s_subgraph.nodes_of_main_graph
        capacities_score = [ get_attribute_node(s_network, s_node, :cap) * 
                sum(get_attribute_edge(s_network, get_edge(s_network, s_node, s_neighbor), :cap) for s_neighbor in neighbors(s_network, s_node)) 
                for s_node in cluster]
        best_node = cluster[findmin(capacities_score)[2]]
        for v_node in v_subgraph.nodes_of_main_graph
            temporary_placement[v_node] = best_node
        end
    end


    # Let's map each sub vn on it's corresponding sub sn
    overall_cost = 0
    node_placement = zeros(Integer, nv(instance.v_network))
    edge_routing = Dict()


    # --- SUBMAPPING TIME!
    for v_subgraph in vn_subgraphs

        s_subgraph = assignment_virtual_substrate_subgraphs[v_subgraph]
        sub_instance = Instance(v_subgraph.graph, s_subgraph.graph)
            

        # Additional costs correspond to the routing of cutting edges on the subproblem.
        # The node placement considered is the temporary placement.
        additional_costs = []
        for v_node in vertices(v_subgraph.graph)
            current_addition_costs = [0 for s_node in vertices(s_subgraph.graph)]
            original_v_node = v_subgraph.nodes_of_main_graph[v_node]

            for v_edge in vn_decompo.v_edges_master
                if src(v_edge) == original_v_node
                    placement_of_dst_node = temporary_placement[dst(v_edge)]
                    for s_node in vertices(s_subgraph.graph)
                        original_s_node = s_subgraph.nodes_of_main_graph[s_node]
                        # We use the basic shortest paths for now, I think it would be quite long and useless to take into account the previous routing
                        current_addition_costs[s_node] += base_shortest_paths.dists[original_s_node, placement_of_dst_node] 
                    end
                end
                if dst(v_edge) == original_v_node
                    placement_of_dst_node = temporary_placement[src(v_edge)]
                    for s_node in vertices(s_subgraph.graph)
                        original_s_node = s_subgraph.nodes_of_main_graph[s_node]
                        # We use the basic shortest paths for now, I think it would be quite long and useless to take into account the previous routing
                        current_addition_costs[s_node] += base_shortest_paths.dists[original_s_node, placement_of_dst_node] 
                    end
                end

            end

            push!(additional_costs, current_addition_costs)
        end

        #println("Additional costs yay! $additional_costs")
        sub_mapping, cost = solve_mepso_custom(sub_instance, additional_costs; nb_particle=25, nb_iter=100, time_max=0.2, print_things=false)
        
        if isnothing(sub_mapping) # invalid submapping!
            println("Submapping failed.")
            return Dict(
                "mapping" => nothing,
                "mapping_cost" => 1e10,
                "solving_time" => time() - time_beginning
            )     
        end
    
        
        real_cost = 0

        for v_node in vertices(v_subgraph.graph)

            original_v_node = v_subgraph.nodes_of_main_graph[v_node]
            original_s_node =  s_subgraph.nodes_of_main_graph[sub_mapping.node_placement[v_node]]
            node_placement[original_v_node] = original_s_node
            temporary_placement[original_v_node] = original_s_node

            real_cost += s_network[original_s_node][:cost]
        end


        for v_edge in edges(v_subgraph.graph)
            original_v_edge = get_edge(v_network, v_subgraph.nodes_of_main_graph[src(v_edge)], v_subgraph.nodes_of_main_graph[dst(v_edge)])
            used_edges = []
            for s_edge in sub_mapping.edge_routing[v_edge].edges
                real_s_edge = get_edge(s_network_dir, s_subgraph.nodes_of_main_graph[src(s_edge)], s_subgraph.nodes_of_main_graph[dst(s_edge)])
                push!(used_edges, real_s_edge)

                real_cost += s_network_dir[src(real_s_edge), dst(real_s_edge)][:cost]
            end
            edge_routing[original_v_edge] = order_path(s_network_dir, used_edges, node_placement[src(original_v_edge)], node_placement[dst(original_v_edge)]) 
        end

        overall_cost += real_cost
    end        


    # --- cutting edge routing time!
    edge_routing, additional_routing_cost = route_cut_edges(instance, vn_decompo, node_placement, edge_routing)
    
    if isnothing(edge_routing) # invalid additional edge routing!
        println("Additional routing failed.")
        return Dict(
            "mapping" => nothing,
            "mapping_cost" => 1e10,
            "solving_time" => time() - time_beginning
        )     
    end


    
    overall_cost += additional_routing_cost


    println("We obtained a mapping of cost $overall_cost !")

    result = Dict()
    result["mapping"] = Mapping(v_network, s_network, node_placement, edge_routing)
    result["mapping_cost"] = overall_cost 
    result["solving_time"] = (time() - time_beginning)
    return result
end







function route_cut_edges(instance, vn_decompo, v_node_placement, edge_routing)
    s_network = instance.s_network

    s_network_dir = instance.s_network_dir
    additional_cost = 0
    s_network_dir_copy = deepcopy(s_network_dir)


    # Capacities and costs on edges
    capacities_edges = Dict{Tuple{Int, Int}, Int}()
    for s_edge in edges(s_network)
        cap = s_network[src(s_edge), dst(s_edge)][:cap]
        capacities_edges[(src(s_edge),dst(s_edge))] = cap
        capacities_edges[(dst(s_edge),src(s_edge))] = cap
    end

    # dstmax matrix
    distmx = zeros(Int, nv(s_network_dir_copy), nv(s_network_dir_copy))
    for s_edge in edges(s_network_dir_copy)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir_copy, s_edge, :cost)
    end

    # remove the edge routing already done from the capacities
    for v_edge in keys(edge_routing)

        path = edge_routing[v_edge]
        
        for s_edge in path.edges
            capacities_edges[(src(s_edge), dst(s_edge))] -= 1
            capacities_edges[(dst(s_edge), src(s_edge))] -= 1
            if capacities_edges[(src(s_edge), dst(s_edge))] <= 0
                rem_edge!(s_network_dir_copy, src(s_edge), dst(s_edge))
                rem_edge!(s_network_dir_copy, dst(s_edge), src(s_edge))
            end
        end

    end
    
    for v_edge in vn_decompo.v_edges_master

        s_src = v_node_placement[src(v_edge)]
        s_dst = v_node_placement[dst(v_edge)]
        shortest_path = a_star(s_network_dir_copy, s_src, s_dst, distmx)

        if shortest_path == []
            #println("No shortest path found: the graph is full!")
            #println("I had the following routing: $edge_routing")
            return nothing, 10e10
        end

        edge_routing[v_edge] = order_path(s_network_dir, shortest_path, s_src, s_dst) 

        for s_edge in shortest_path
            capacities_edges[(src(s_edge), dst(s_edge))] -= 1
            capacities_edges[(dst(s_edge), src(s_edge))] -= 1
            if capacities_edges[(src(s_edge), dst(s_edge))] <= 0
                rem_edge!(s_network_dir_copy, src(s_edge), dst(s_edge))
                rem_edge!(s_network_dir_copy, dst(s_edge), src(s_edge))
            end
        end

        additional_cost += edge_routing[v_edge].cost
    end


    return edge_routing, additional_cost

end






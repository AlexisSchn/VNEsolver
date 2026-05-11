using Random
using Graphs, MetaGraphsNext


function solve_local_search_pricer(sub_instance, dual_costs; nb_iterations=500)



    function shortest_path_routing(v_node_placement, get_routing=false)

        capacities_edges_copy = copy(capacities_edges)
        s_network_copy_dir_copy_ofgraph = nothing
        is_still_original_s_network = true

        edge_routing = Dict{Edge, Path}()
        overall_routing_costs = 0

        for v_edge in edges(v_network)

            s_src = v_node_placement[src(v_edge)]
            s_dst = v_node_placement[dst(v_edge)]

            nodes_of_path = Int[]
            cost_of_routing_current_edge = 0

            # Get the shortest path
            if is_still_original_s_network
                nodes_of_path = [s_dst]
                v = s_dst
                while v != s_src
                    u = shortest_paths.parents[s_src, v]
                    push!(nodes_of_path, u)
                    v = shortest_paths.parents[s_src, v]
                end
                reverse!(nodes_of_path)
                cost_of_routing_current_edge = shortest_paths.dists[s_src, s_dst]
            else # A smarter thing to do here, would be to use the basic path, and if it's using a removed edge, then compute the astar path. But it does not cost so much time so.
                edges_of_path = a_star(s_network_copy_dir_copy_ofgraph, s_src, s_dst, distmx)
                if edges_of_path == [] # No paths found: the graph is full!
                    return Dict(), 10e9
                end
                push!(nodes_of_path, src(edges_of_path[1]))
                for edge in edges_of_path
                    push!(nodes_of_path, dst(edge))
                    cost_of_routing_current_edge += distmx[src(edge), dst(edge)]
                end
            end
            
            overall_routing_costs += cost_of_routing_current_edge

            # Removing the capacities. If no more capacities, removing the edge.
            for i_node in 1:length(nodes_of_path)-1
                current_src = nodes_of_path[i_node]
                current_dst = nodes_of_path[i_node+1]
                capacities_edges_copy[(current_src), (current_dst)] -= 1
                capacities_edges_copy[(current_dst), (current_src)] -= 1

                if capacities_edges_copy[(current_src), (current_dst)] == 0
                    if is_still_original_s_network
                        s_network_copy_dir_copy_ofgraph = deepcopy(s_network_dir.graph)
                        is_still_original_s_network = false
                    end
                    rem_edge!(s_network_copy_dir_copy_ofgraph, current_src, current_dst)
                    rem_edge!(s_network_copy_dir_copy_ofgraph, current_dst, current_src)
                end
            end

            # Sometime it is necessary to get the actual routing, but not all the time. 
            # Since it takes some time to retrieve the edges, due to poor design choices, I only do it when necessary.
            if get_routing
                edges_of_path = Edge[]
                real_cost_of_path = 0

                for i_node in 1:length(nodes_of_path)-1
                    current_src = nodes_of_path[i_node]
                    current_dst = nodes_of_path[i_node+1]
                    push!(edges_of_path, get_edge(s_network_dir, current_src, current_dst))
                    real_cost_of_path += s_network_dir[current_src, current_dst][:cost]
                end
                path = Path(s_src, s_dst, edges_of_path, real_cost_of_path)
                edge_routing[v_edge] = path
            end
        end

        return edge_routing, overall_routing_costs
    end



    function complete_partial_placement(partial_placement)

        already_placed_v_nodes = Int[]
        for v_node in 1:nv(v_network)
            if partial_placement[v_node] != 0
                push!(already_placed_v_nodes, v_node)
            end
        end

        placement = copy(partial_placement)

        # Extreme case where no nodes are kept: we start from a random node again
        if isempty(already_placed_v_nodes)
            keep_on = true
            while keep_on 
                s_node = rand(1:nv(s_network))
                v_node = rand(1:nv(v_network))
                if node_capacities[s_node] >= 1
                    placement[v_node] = s_node
                    push!(already_placed_v_nodes, v_node)
                    keep_on = false
                end
            end
        end

        next_v_nodes = Int[]
        for v_node in already_placed_v_nodes
            for v_neighbor in neighbors(v_network, v_node)
                if v_neighbor ∉ next_v_nodes && v_neighbor ∉ already_placed_v_nodes
                    push!(next_v_nodes, v_neighbor)
                end
            end
        end

        while !isempty(next_v_nodes)

            # Take a node of the list
            shuffle!(next_v_nodes)
            v_node = popfirst!(next_v_nodes)

            # Get neighbors already placed
            placement_neighbors = Int[]
            for v_neighbor in neighbors(v_network, v_node)
                if v_neighbor ∈ already_placed_v_nodes
                    push!(placement_neighbors, placement[v_neighbor])
                end
            end

            # Choose 10 random nodes (or less if the network is not that big)
            number_s_nodes = min(floor(nv(s_network)/5), 10)
            some_s_nodes = Int[]
            while length(some_s_nodes) < number_s_nodes
                s_node = rand(1:nv(s_network))
                if (node_capacities[s_node] >= 1) && (s_node ∉ placement)
                    push!(some_s_nodes, s_node)
                end
            end

            # Rank them according distance to already placed nodes and capacity
            distances = [ sum(shortest_paths.dists[s_src, s_node] for s_src in placement_neighbors) for s_node in some_s_nodes]
            distances_norm = (distances .- minimum(distances)) ./ (maximum(distances) - minimum(distances) + 1e-9)

            capacities = [ s_node_scores[s_node] for s_node in some_s_nodes]
            capacities_norm = (capacities .- minimum(capacities)) ./ (maximum(capacities) - minimum(capacities) + 1e-9)

            costs = [node_costs[s_node] for s_node in some_s_nodes]
            costs_norm = (costs .- minimum(costs)) ./ (maximum(costs) - minimum(costs) + 1e-9)

            final_scores = 1. .* distances_norm .+ 1. .* costs_norm .+ 0.5 .* ( 1 .-capacities_norm)

            selected_idx = argmin(final_scores)
            s_node_selected = some_s_nodes[selected_idx]


            # Finish the work
            placement[v_node] = s_node_selected
            push!(already_placed_v_nodes, v_node)
            for v_neighbor in neighbors(v_network, v_node)
                if v_neighbor ∉ already_placed_v_nodes && v_neighbor ∉ next_v_nodes 
                    push!(next_v_nodes, v_neighbor)
                end
            end

        end

        placement_cost = 0
        for v_node in vertices(v_network)
            placement_cost += node_costs[placement[v_node]]
        end

        return placement, placement_cost


    end



    v_network = sub_instance.v_network
    s_network = sub_instance.s_network
    s_network_dir = sub_instance.s_network_dir


    bad_result = (sub_mapping=nothing,
        real_cost =10e9,
        reduced_cost=10e9
    )


    #---- Make sure there are enough capacited nodes
    nodes_with_caps = 0
    for s_node in vertices(s_network)
        if s_network[s_node][:cap] >= 1
            nodes_with_caps += 1
        end
    end
    if nodes_with_caps < nv(v_network)
        #println("What the hell? Not enough capacited nodes...")
        return bad_result
    end

        
    #---- Usefull things
    node_capacities = [get_attribute_node(s_network, s_node, :cap) for s_node in vertices(s_network)]
    node_costs = [get_attribute_node(s_network, s_node, :cost) - dual_costs.capacity_s_node[s_node] for s_node in vertices(s_network)]
    capacities_edges = Dict{Tuple{Int, Int}, Int}()
    for s_edge in edges(s_network)
        cap = s_network[src(s_edge), dst(s_edge)][:cap]
        capacities_edges[(src(s_edge),dst(s_edge))] = cap
        capacities_edges[(dst(s_edge),src(s_edge))] = cap
    end

    # Greedy score based on capacities for nodes
    s_node_scores = [ node_capacities[s_node]  * 
                        sum(capacities_edges[(s_node,s_neighbor)] for s_neighbor in neighbors(s_network, s_node)) 
                    for s_node in vertices(s_network)
    ]



    # shortest path of the substrate network --- WITH DUAL COSTS!
    distmx = zeros(nv(s_network), nv(s_network))
    for s_edge in edges(s_network)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network, s_edge, :cost) - dual_costs.capacity_s_edge[s_edge]
        distmx[dst(s_edge), src(s_edge)] = get_attribute_edge(s_network, s_edge, :cost) - dual_costs.capacity_s_edge[s_edge]
    end
    shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)

    


    # Getting the most central virtual node
    most_central_v_node = argmin(closeness_centrality(v_network))

    # Loops related things
    best_cost = 10e9
    best_placement = zeros(nv(v_network))
    
    

    # Construct initial mapping
    # We iterate, for the case where it's hard to find a feasible one

    iter = 0
    do_init = true
    while do_init && iter < nb_iterations

        s_node = rand(1:nv(s_network))

        placement = zeros(Int32, nv(v_network))
        placement[most_central_v_node] = s_node
    
        placement, placement_cost = complete_partial_placement(placement) 
        routing, routing_cost = shortest_path_routing(placement)
    
        best_cost = placement_cost + routing_cost
        best_placement = placement

        if best_cost < 10e6
            do_init = false
        end

        iter += 1
    end


    # Local search iterations

    while iter < nb_iterations

        placement = copy(best_placement)
        
        # delete some random star in the virtual network
        some_v_node = rand(1:nv(v_network))
        v_nodes_deleted = [some_v_node]
        for v_neighbor in neighbors(v_network, some_v_node) 
            push!(v_nodes_deleted, v_neighbor)
        end

        for v_node in v_nodes_deleted
            placement[v_node] = 0
        end
            

        # reconstruct
        placement, placement_cost = complete_partial_placement(placement) 
        routing, routing_cost = shortest_path_routing(placement)


        current_cost = placement_cost + routing_cost
        if current_cost < best_cost
            best_cost = current_cost
            best_placement = placement
        end

        iter += 1
    end



    if best_cost>10e6
        return bad_result
    end

    routing, routing_cost= shortest_path_routing(best_placement, true)
    best_mapping = Mapping(v_network, s_network, best_placement, routing)

    real_cost = get_cost_placement(best_mapping) + get_cost_routing(best_mapping)

    return (sub_mapping=best_mapping,
            real_cost=real_cost,
            reduced_cost=best_cost
    )

end
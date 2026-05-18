using Random
using Graphs, MetaGraphsNext
using StatsBase

includet("../utils/import_utils.jl")




function solve_greedy(instance; nb_greedy = 10, time_max = 10)



    
    function shortest_path_routing(v_node_placement; get_routing=false)

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
                        s_network_copy_dir_copy_ofgraph = deepcopy(instance.s_network_dir.graph)
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
                for i_node in 1:length(nodes_of_path)-1
                    current_src = nodes_of_path[i_node]
                    current_dst = nodes_of_path[i_node+1]
                    push!(edges_of_path, get_edge(s_network_dir, current_src, current_dst))
                end
                path = Path(s_src, s_dst, edges_of_path, cost_of_routing_current_edge)
                edge_routing[v_edge] = path
            end
        end

        return edge_routing, overall_routing_costs
    end


    function shortest_path_routing_bruh(v_node_placement; get_routing=false)

        capacities_edges_copy = copy(capacities_edges)
        s_network_copy_dir_copy_ofgraph = deepcopy(instance.s_network_dir.graph)

        edge_routing = Dict{Edge, Path}()
        overall_routing_costs = 0

        for v_edge in edges(v_network)

            s_src = v_node_placement[src(v_edge)]
            s_dst = v_node_placement[dst(v_edge)]
            demand_curr_edge = edge_demands[(src(v_edge), dst(v_edge))]

            nodes_of_path = Int[]
            cost_of_routing_current_edge = 0

            # Routing
            #!!!!!!! TODO CHANGEEEEEEE IT FOR NON-UNIFORM DEMANDS
            edges_of_path = a_star(s_network_copy_dir_copy_ofgraph, s_src, s_dst, distmx)
            if edges_of_path == [] # No paths found: the graph is full!
                return Dict(), 10e9
            end
            push!(nodes_of_path, src(edges_of_path[1]))
            for edge in edges_of_path
                push!(nodes_of_path, dst(edge))
                cost_of_routing_current_edge += distmx[src(edge), dst(edge)] * demand_curr_edge
            end
            
            overall_routing_costs += cost_of_routing_current_edge

            # Removing the capacities. If no more capacities, removing the edge.
            for i_node in 1:length(nodes_of_path)-1
                current_src = nodes_of_path[i_node]
                current_dst = nodes_of_path[i_node+1]
                capacities_edges_copy[(current_src), (current_dst)] -= demand_curr_edge
                capacities_edges_copy[(current_dst), (current_src)] -= demand_curr_edge

                if capacities_edges_copy[(current_src), (current_dst)] == 0
                    rem_edge!(s_network_copy_dir_copy_ofgraph, current_src, current_dst)
                    rem_edge!(s_network_copy_dir_copy_ofgraph, current_dst, current_src)
                end
            end

            # Sometime it is necessary to get the actual routing, but not all the time. 
            # Since it takes some time to retrieve the edges, due to poor design choices, I only do it when necessary.
            if get_routing
                edges_of_path = Edge[]
                for i_node in 1:length(nodes_of_path)-1
                    current_src = nodes_of_path[i_node]
                    current_dst = nodes_of_path[i_node+1]
                    push!(edges_of_path, get_edge(s_network_dir, current_src, current_dst))
                end
                path = Path(s_src, s_dst, edges_of_path, cost_of_routing_current_edge)
                edge_routing[v_edge] = path
            end
        end

        return edge_routing, overall_routing_costs
    end




    function complete_partial_placement(partial_placement; nb_nodes_to_try=10)

        placement = copy(partial_placement)

        #println("Partial placement : $partial_placement")
        already_placed_v_nodes = filter(v_node -> placement[v_node] != 0, vertices(v_network))

        # Extreme case where no nodes are kept: we start from a random node again
        if isempty(already_placed_v_nodes)
            placement[most_central_v_node] = rand(capacited_nodes)
            push!(already_placed_v_nodes, most_central_v_node)
        end

        next_v_nodes = setdiff(
            union([neighbors(v_network, v) for v in already_placed_v_nodes]...),
            already_placed_v_nodes,
        )

        possible_s_nodes = filter(s_node -> s_node ∉ placement, capacited_nodes)

        while !isempty(next_v_nodes)

            
            # Take a node of the list
            shuffle!(next_v_nodes)
            v_node = popfirst!(next_v_nodes)

            # Get neighbors already placed
            placement_neighbors = [placement[s_neigh] for s_neigh in filter(v_neighbor -> placement[v_neighbor] != 0, neighbors(v_network, v_node))]

            # Choose 10 random nodes (or less if the network is not that big)
            number_s_nodes = min(length(possible_s_nodes), nb_nodes_to_try)

            if number_s_nodes == 0
                return [], 10e9
            end

            some_s_nodes = sample(possible_s_nodes, number_s_nodes; replace=false)


            # NODE RANKING AND SELECTION
            distances = [ sum(shortest_paths.dists[s_src, s_node] for s_src in placement_neighbors) for s_node in some_s_nodes]
            distances_norm = (distances .- minimum(distances)) ./ (maximum(distances) - minimum(distances) + 1e-9)

            capacities = [ s_node_scores[s_node] for s_node in some_s_nodes]
            capacities_norm = (capacities .- minimum(capacities)) ./ (maximum(capacities) - minimum(capacities) + 1e-9)

            costs = [ node_costs[s_node] for s_node in some_s_nodes]
            costs_norm = (costs .- minimum(costs)) ./ (maximum(costs) - minimum(costs) + 1e-9)

            final_scores = 5. .* distances_norm .+ 0 .* costs_norm .+ 0. .* ( 1 .-capacities_norm)
            selected_idx = argmin(final_scores)
            s_node_selected = some_s_nodes[selected_idx]


            # Finish the work
            placement[v_node] = s_node_selected
            push!(already_placed_v_nodes, v_node)
            next_v_nodes = next_v_nodes ∪ filter(v_neighbor->v_neighbor ∉ already_placed_v_nodes, neighbors(v_network, v_node) )
            possible_s_nodes = filter(!=(s_node_selected), possible_s_nodes)
        end

        placement_cost = 0
        for v_node in vertices(v_network)
            placement_cost += node_costs[placement[v_node]] 
        end

        return placement, placement_cost

    end

    function complete_partial_placement_bruh(partial_placement; nb_nodes_to_try=10)

        placement = copy(partial_placement)

        #println("Partial placement : $partial_placement")
        already_placed_v_nodes = filter(v_node -> placement[v_node] != 0, vertices(v_network))


        next_v_nodes = setdiff(
            union([neighbors(v_network, v) for v in already_placed_v_nodes]...),
            already_placed_v_nodes,
        )

        possible_s_nodes = filter(s_node -> s_node ∉ placement, vertices(s_network))

        while !isempty(next_v_nodes)

            # Take a node of the list
            shuffle!(next_v_nodes)
            v_node = popfirst!(next_v_nodes)

            # Get neighbors already placed
            placement_neighbors = [placement[s_neigh] for s_neigh in filter(v_neighbor -> placement[v_neighbor] != 0, neighbors(v_network, v_node))]

            # Choose some random nodes
            number_s_nodes = min(length(possible_s_nodes), nb_nodes_to_try)
            some_s_nodes = sample(possible_s_nodes, number_s_nodes; replace=false)


            # NODE RANKING AND SELECTION
            distances = [ sum(shortest_paths.dists[s_src, s_node] for s_src in placement_neighbors) for s_node in some_s_nodes]
            distances_norm = (distances .- minimum(distances)) ./ (maximum(distances) - minimum(distances) + 1e-9)

            #capacities = [ s_node_scores[s_node] for s_node in some_s_nodes]
            #capacities_norm = (capacities .- minimum(capacities)) ./ (maximum(capacities) - minimum(capacities) + 1e-9)

            #costs = [ node_costs[s_node] for s_node in some_s_nodes]
            #costs_norm = (costs .- minimum(costs)) ./ (maximum(costs) - minimum(costs) + 1e-9)

            final_scores = [ node_capacities[some_s_nodes[idx_node]] >= node_demands[v_node] ? distances_norm[idx_node] : 10e9
                                for idx_node in 1:length(some_s_nodes)]
            selected_idx = argmin(final_scores)
            s_node_selected = some_s_nodes[selected_idx]

            if minimum(final_scores) > 10e8
                return [], 10e9
            end

            # Finish the work
            placement[v_node] = s_node_selected
            push!(already_placed_v_nodes, v_node)
            next_v_nodes = next_v_nodes ∪ filter(v_neighbor->v_neighbor ∉ already_placed_v_nodes, neighbors(v_network, v_node) )
            possible_s_nodes = filter(!=(s_node_selected), possible_s_nodes)
        end

        placement_cost = 0
        for v_node in vertices(v_network)
            placement_cost += node_costs[placement[v_node]] * node_demands[v_node]
        end

        return placement, placement_cost

    end



    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir




    #---- Usefull things
    node_demands = [get_attribute_node(v_network, v_node, :dem) for v_node in vertices(v_network)]
    edge_demands = Dict{Tuple{Int, Int}, Int}()
    for v_edge in edges(v_network)
        edge_demands[(src(v_edge),dst(v_edge))] = get_attribute_edge(v_network, v_edge, :dem)
    end
    

    node_capacities = [get_attribute_node(s_network, s_node, :cap) for s_node in vertices(s_network)]
    node_costs = [get_attribute_node(s_network, s_node, :cost) for s_node in vertices(s_network)]
    capacities_edges = Dict{Tuple{Int, Int}, Int}()
    for s_edge in edges(s_network)
        cap = s_network[src(s_edge), dst(s_edge)][:cap]
        capacities_edges[(src(s_edge),dst(s_edge))] = cap
        capacities_edges[(dst(s_edge),src(s_edge))] = cap
    end
    capacited_nodes = [s_node for s_node in vertices(s_network) if node_capacities[s_node] ≥ 1]

    
    # Greedy score based on capacities for nodes
    s_node_scores = [ sum(capacities_edges[(s_node,s_neighbor)] + node_capacities[s_neighbor] for s_neighbor in neighbors(s_network, s_node) ) 
                    for s_node in vertices(s_network)]


    # shortest path of the substrate network
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network_dir)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir, s_edge, :cost)
    end
    shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)



    time_beginning = time()
    best_cost = 10e8
    best_mapping = nothing
    time_overall = 0
    iter = 1
    centrality_nodes = closeness_centrality(s_network)
    most_central_v_node = argmin(closeness_centrality(v_network))

    while iter <= nb_greedy && time()-time_beginning<time_max

        # Getting the most central virtual node
            
        # Construct initial mapping
        init_s_nodes_scores = [ node_capacities[s_node] >= node_demands[most_central_v_node] ? centrality_nodes[s_node] + 1. * rand() : 10e9
                                    for s_node in vertices(s_network) ]
        
        s_node_start = argmin(init_s_nodes_scores)
        if minimum(init_s_nodes_scores) > 10e8
            return (mapping=nothing,
                    mapping_cost=10e9
            )
        end


        placement = zeros(Int32, nv(v_network))
        placement[most_central_v_node] = s_node_start

        
        placement, placement_cost = complete_partial_placement(placement; nb_nodes_to_try=nv(s_network)) 
        if placement_cost < 10e8
            routing, routing_cost = shortest_path_routing(placement, get_routing=false)
        else
            routing_cost = 10e9
        end

        total_cost = placement_cost + routing_cost

        if total_cost < best_cost
            best_cost = total_cost
            routing, routing_cost = shortest_path_routing(placement, get_routing=true)
            best_mapping = Mapping(v_network, s_network, placement, routing)
            println("New best mapping with cost $best_cost")
        end

        iter += 1
    end

    println("Found $best_cost with $iter iterations and $(time()-time_beginning) time")
    return (mapping=best_mapping, mapping_cost=best_cost)
end



using Random
using Graphs, MetaGraphsNext
using StatsBase

includet("../../utils/import_utils.jl")

function solve_greedy(instance, given_placement)



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



    function complete_partial_placement(partial_placement; nb_nodes_to_try=10)

        placement = copy(partial_placement)

        #println("Partial placement : $partial_placement")
        already_placed_v_nodes = filter(v_node -> placement[v_node] != 0, vertices(v_network))


        next_v_nodes = setdiff(
            union([neighbors(v_network, v) for v in already_placed_v_nodes]...),
            already_placed_v_nodes,
        )

        possible_s_nodes = filter(s_node -> s_node ∉ placement, capacited_nodes)
        #println("next v nodes: $next_v_nodes")
        while !isempty(next_v_nodes)
            #println("hum")
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

            capacities = [ s_node_scores[s_node] for s_node in some_s_nodes]
            capacities_norm = (capacities .- minimum(capacities)) ./ (maximum(capacities) - minimum(capacities) + 1e-9)

            costs = [ node_costs[s_node] for s_node in some_s_nodes]
            costs_norm = (costs .- minimum(costs)) ./ (maximum(costs) - minimum(costs) + 1e-9)

            final_scores = 3. .* distances_norm .+ 1. .* costs_norm .+ 0.5 .* ( 1 .-capacities_norm)
            selected_idx = argmin(final_scores)
            s_node_selected = some_s_nodes[selected_idx]


            # Finish the work
            placement[v_node] = s_node_selected
            push!(already_placed_v_nodes, v_node)
            next_v_nodes = next_v_nodes ∪ filter(v_neighbor->v_neighbor ∉ already_placed_v_nodes, neighbors(v_network, v_node) )
            possible_s_nodes = filter(!=(s_node_selected), possible_s_nodes)
        end

        #println("Yo placement: $placement")
        placement_cost = 0
        for v_node in vertices(v_network)
            placement_cost += node_costs[placement[v_node]] 
        end

        return placement, placement_cost

    end



    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir



        
    #---- Usefull things
    node_capacities = [get_attribute_node(s_network, s_node, :cap) for s_node in vertices(s_network)]
    node_costs = [get_attribute_node(s_network, s_node, :cost) for s_node in vertices(s_network)]
    capacities_edges = Dict{Tuple{Int, Int}, Int}()
    for s_edge in edges(s_network)
        cap = s_network[src(s_edge), dst(s_edge)][:cap]
        capacities_edges[(src(s_edge),dst(s_edge))] = cap
        capacities_edges[(dst(s_edge),src(s_edge))] = cap
    end




    #---- Make sure there are enough capacited nodes
    capacited_nodes = [s_node for s_node in vertices(s_network) if node_capacities[s_node] ≥ 1]

    if length(capacited_nodes) < nv(v_network)
        println("What the hell? Not enough capacited nodes...")
        return  (mapping=nothing, 
                mapping_cost=10e9
        )
    end
    
    # Greedy score based on capacities for nodes
    s_node_scores = [ sum(capacities_edges[(s_node,s_neighbor)] + node_capacities[s_neighbor] for s_neighbor in neighbors(s_network, s_node) ) 
                    for s_node in vertices(s_network)]


    # shortest path of the substrate network
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network_dir)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir, s_edge, :cost)
    end
    shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)


    # Getting the most central virtual node
    most_central_v_node = argmin(closeness_centrality(v_network))
        
    # Construct initial mapping
    placement = zeros(Int32, nv(v_network))
    if isempty(given_placement)
        centrality_nodes = closeness_centrality(s_network)
        s_nodes_scores = [ (centrality_nodes[s_node] + 0.5 * rand() ) for s_node in capacited_nodes ]
        s_node_start = capacited_nodes[argmin(s_nodes_scores)]
        placement[most_central_v_node] = s_node_start
    else
        for v_node in keys(given_placement)
            placement[v_node] = given_placement[v_node]
        end
    end
    placement, placement_cost = complete_partial_placement(placement; nb_nodes_to_try=length(capacited_nodes)) 

    routing, routing_cost = shortest_path_routing(placement, true)
    total_cost = placement_cost + routing_cost

    if total_cost>10e6
        return (mapping=nothing,
                    mapping_cost=10e9
        )
    end

    final_mapping = Mapping(v_network, s_network, placement, routing)

    #println("Find the solution of $total_cost in $(time()-time_beginning)")
    return (mapping = final_mapping,
            mapping_cost = total_cost)
end


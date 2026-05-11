

using Graphs, MetaGraphsNext


function solve_mepso_custom(instance, additional_costs; nb_particle=25, nb_iter=50, time_max=5, print_things=true)



    function shortest_path_routing(v_node_placement, get_routing=false)

        capacities_edges_copy = deepcopy(capacities_edges)
        s_network_copy_dir_copy_ofgraph = nothing
        is_still_original_s_network = true

        edge_routing = Dict()
        overall_routing_costs = 0

        for v_edge in edges(v_network)

            s_src = v_node_placement[src(v_edge)]
            s_dst = v_node_placement[dst(v_edge)]

            nodes_of_path = []
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
                edges_of_path = []
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



    function construct_greedy_random_placement()

        placement = zeros(Int, nv(v_network))

        next_v_nodes = []
        already_placed_v_nodes = []

        keep_on = true
        while keep_on 
            s_node = rand(1:nv(s_network))
            v_node = rand(1:nv(v_network))
            if node_capacities[s_node] >= 1
                placement[v_node] = s_node
                push!(already_placed_v_nodes, v_node)
                append!(next_v_nodes, neighbors(v_network, v_node))
                keep_on = false
            end
        end

        while !isempty(next_v_nodes)

            # Take a node of the list
            v_node = popfirst!(next_v_nodes)
            # Get neighbors already placed
            placement_neighbors = []
            for v_neighbor in neighbors(v_network, v_node)
                if v_neighbor ∈ already_placed_v_nodes
                    push!(placement_neighbors, placement[v_neighbor])
                end
            end

            # Choose 10 random nodes (or less if the network is not that big)
            number_s_nodes = min(floor(nv(s_network)/5), 10)
            some_s_nodes = []
            while length(some_s_nodes) < 7
                s_node = rand(1:nv(s_network))
                if (node_capacities[s_node] >= 1) && (s_node ∉ placement)
                    push!(some_s_nodes, s_node)
                end
            end

            # Rank them according distance to already placed nodes and capacity
            distances = [ (sum(shortest_paths.dists[s_src, s_node] for s_src in placement_neighbors)
                            + additional_costs[v_node][s_node])
                            for s_node in some_s_nodes]
            distances_norm = (distances .- minimum(distances)) ./ (maximum(distances) - minimum(distances) + 1e-9)

            capacities = [ s_node_scores[s_node] for s_node in some_s_nodes]
            capacities_norm = (capacities .- minimum(capacities)) ./ (maximum(capacities) - minimum(capacities) + 1e-9)

            final_scores = 1. .* distances_norm .+ 1. .* ( 1 .-capacities_norm) + 0.5 * rand(length(some_s_nodes))

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
            placement_cost += (node_costs[placement[v_node]] + additional_costs[v_node][placement[v_node]])
        end

        routing, routing_cost = shortest_path_routing(placement)

        overall_cost = placement_cost + routing_cost

        return placement, overall_cost

    end



    function minus(pos1, pos2)

        res=[]
        for i in 1:length(pos1)
            if pos1[i] == pos2[i]
                push!(res, 1)
            else
                push!(res, 0)
            end
        end
        return res
    end
    
    
    function plus(vel_inertia, vel_pb, vel_gb)
    
        p_inertia = 0.2
        p_attraction_personal = 0.3
        p_attraction_global = 0.5
    
        new_velocity = []
        for i in 1:length(vel_inertia)
            r = rand()
            if r < p_inertia
                push!(new_velocity, vel_inertia[i])
            elseif r < (p_inertia + p_attraction_personal)
                push!(new_velocity, vel_pb[i])
            else
                push!(new_velocity, vel_gb[i])
            end
        end
    
        return new_velocity

    end

    
    function times_greedy_random(position, velocity)


        # Get the partial placement (destruction from PSO)
        placement = zeros(Int, nv(v_network))

        already_placed_v_nodes = []
        for v_node in 1:nv(instance.v_network)
            if velocity[v_node] == 1
                placement[v_node] = position[v_node]
                push!(already_placed_v_nodes, v_node)
            end
        end

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

        next_v_nodes = []
        for v_node in already_placed_v_nodes
            for v_neighbor in neighbors(v_network, v_node)
                if v_neighbor ∉ next_v_nodes && v_neighbor ∉ already_placed_v_nodes
                    push!(next_v_nodes, v_neighbor)
                end
            end
        end

        while !isempty(next_v_nodes)

            # Take a node of the list
            v_node = popfirst!(next_v_nodes)
            # Get neighbors already placed
            placement_neighbors = []
            for v_neighbor in neighbors(v_network, v_node)
                if v_neighbor ∈ already_placed_v_nodes
                    push!(placement_neighbors, placement[v_neighbor])
                end
            end

            # Choose 10 random nodes (or less if the network is not that big)
            number_s_nodes = min(floor(nv(s_network)/5), 10)
            some_s_nodes = []
            while length(some_s_nodes) < 7
                s_node = rand(1:nv(s_network))
                if (node_capacities[s_node] >= 1) && (s_node ∉ placement)
                    push!(some_s_nodes, s_node)
                end
            end

            # Rank them according distance to already placed nodes and capacity
            distances = [ (sum(shortest_paths.dists[s_src, s_node] for s_src in placement_neighbors)
                            + additional_costs[v_node][s_node])
                            for s_node in some_s_nodes]
            distances_norm = (distances .- minimum(distances)) ./ (maximum(distances) - minimum(distances) + 1e-9)

            capacities = [ s_node_scores[s_node] for s_node in some_s_nodes]
            capacities_norm = (capacities .- minimum(capacities)) ./ (maximum(capacities) - minimum(capacities) + 1e-9)

            final_scores = 1. .* distances_norm .+ 1. .* ( 1 .-capacities_norm) + 0.5 * rand(length(some_s_nodes))

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
            placement_cost += (node_costs[placement[v_node]] + additional_costs[v_node][placement[v_node]])
        end

        return placement, placement_cost
    end


    time_start = time()
    time_sp = 0
    time_pos = 0

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    #---- Make sure there are enough capacited nodes
    nodes_with_caps = 0
    for s_node in vertices(s_network)
        if s_network[s_node][:cap] >= 1
            nodes_with_caps += 1
        end
    end
    if nodes_with_caps < nv(v_network)
        return nothing, 10e9
    end



    #---- usefull things for the resolution
    node_capacities = [get_attribute_node(s_network, s_node, :cap) for s_node in vertices(s_network)]
    node_costs = [get_attribute_node(s_network, s_node, :cost) for s_node in vertices(s_network)]
    capacities_edges = Dict{Tuple{Int, Int}, Int}()
    for s_edge in edges(s_network)
        cap = s_network[src(s_edge), dst(s_edge)][:cap]
        capacities_edges[(src(s_edge),dst(s_edge))] = cap
        capacities_edges[(dst(s_edge),src(s_edge))] = cap
    end

    # greedy stuff
    s_node_scores = [ get_attribute_node(s_network, s_node, :cap) * sum(get_attribute_edge(s_network, get_edge(s_network, s_node, s_neighbor), :cap) for s_neighbor in neighbors(s_network, s_node)) 
            for s_node in vertices(s_network)]


    # shortest path things
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network_dir)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir, s_edge, :cost)
    end
    shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)

    # ------ things for the PSO algorithm
    position = []
    velocity = []

    personal_best = []
    personal_best_cost = []

    global_best = nothing
    global_best_cost = 10e9

    # initialization
    print_things && print("initialization... ")

    for particle in 1:nb_particle

        placement, overall_cost = construct_greedy_random_placement()

        push!(position, placement)
        push!(personal_best, placement)
        push!(personal_best_cost, overall_cost)

        if overall_cost < global_best_cost
            global_best = position[particle]
            global_best_cost = overall_cost
            print_things && println("We got a new best solution! value $overall_cost")
        end

        push!(velocity, ones(nv(v_network)))
    end
    print_things && println("Initialization done, Starting iterations...")

    # iterations
    iter = 1
    time_total = 0
    while iter < nb_iter && time_total < time_max
        for particle in 1:nb_particle

            if personal_best_cost[particle] > 99999 # if the first isnt good, we reinitialized
                
                placement, overall_cost = construct_greedy_random_placement()
        
                
                if overall_cost < 999999
                    position[particle] = placement
                    personal_best[particle] = placement
                    personal_best_cost[particle] = overall_cost
                end
        
                if overall_cost < global_best_cost
                    global_best = position[particle]
                    global_best_cost = overall_cost
                    print_things && println("We got a new best solution! value $overall_cost")
                end



            else # we do a normal iteration
                velocity[particle] = plus( velocity[particle], 
                                            minus(personal_best[particle], position[particle]), 
                                            minus(global_best, position[particle]))
                position[particle], placement_cost = times_greedy_random(position[particle], velocity[particle])

                time_beg_sp = time()
                routing, routing_cost = shortest_path_routing(position[particle])
                time_sp += time()-time_beg_sp
                overall_cost = placement_cost + routing_cost
            end
            
            if overall_cost < personal_best_cost[particle]
                    personal_best[particle] = position[particle]
                    personal_best_cost[particle] = overall_cost
            end
            if overall_cost < global_best_cost
                    global_best_cost = overall_cost
                    global_best = position[particle]
                    print_things && println("We got a new best solution! value $global_best_cost")
            end

        end

        iter += 1
        time_total = time() - time_start

    end

    #println("Final best solution: $global_best")
    print_things && println("PSO finished at iteration $nb_iter, finished in $(time()-time_start)s, best solution: $global_best_cost")

    if global_best_cost > 10e8
        return nothing, 10e9
    end

    routing, routing_cost= shortest_path_routing(global_best, true)
    final_mapping = Mapping(v_network, s_network, global_best, routing)

    return final_mapping, global_best_cost
end



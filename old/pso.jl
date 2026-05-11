
using Graphs, MetaGraphsNext
includet("../utils/import_utils.jl")




function solve_PSO(instance; nb_particle=25, nb_iter=50, time_max=5, print_things=true)



    function shortest_path_routing(v_node_placement, get_routing=false) # TODO: do for real the thing where you remove the edges bro. It's important......

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



    
    function construct_random_placement()

        placement = []
        placement_cost=0
        for i in 1:nv(v_network)
            keep_on = true
            while keep_on                
                s_node = rand(1:nv(s_network))
                if s_node ∉ placement && node_capacities[s_node]>0
                    push!(placement, s_node)
                    keep_on=false
                    placement_cost+= node_costs[s_node]
                end
            end
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

    

    function times(position, velocity)

        new_placement = []
        placement_cost = 0
    
        for i in 1:nv(instance.v_network)
            if velocity[i] == 1
                push!(new_placement, position[i])
                placement_cost += node_costs[position[i]]
            else
                push!(new_placement, -1)
            end
        end
    
        for i in 1:nv(instance.v_network)
            if new_placement[i] == -1
                keep_on = true
                while keep_on
                    s_node = rand(1:nv(s_network))
                    if s_node ∉ new_placement && node_capacities[s_node]>0
                        new_placement[i]=s_node
                        keep_on=false
                        placement_cost+= node_costs[s_node]
                    end
                end
            end
        end
    
        return new_placement, placement_cost
    end
    

    time_start = time()
    time_sp = 0
    time_pos = 0

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    #---- usefull things for the resolution
    node_capacities = [get_attribute_node(s_network, s_node, :cap) for s_node in vertices(s_network)]
    node_costs = [get_attribute_node(s_network, s_node, :cost) for s_node in vertices(s_network)]
    capacities_edges = Dict{Tuple{Int, Int}, Int}()
    for s_edge in edges(s_network)
        cap = s_network[src(s_edge), dst(s_edge)][:cap]
        capacities_edges[(src(s_edge),dst(s_edge))] = cap
        capacities_edges[(dst(s_edge),src(s_edge))] = cap
    end

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
    global_best_cost = 9999999

    # initialization
    print_things && print("initialization... ")

    for particle in 1:nb_particle

        placement, overall_cost = construct_random_placement()

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
                
                placement, overall_cost = construct_random_placement()
        
                
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
                position[particle], placement_cost = times(position[particle], velocity[particle])

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
    print_things && println("We spent $time_sp s in shortest paths")
    routing, routing_cost_shortest_path = shortest_path_routing(global_best, true)
    final_mapping = Mapping(v_network, s_network, global_best, routing)

    return final_mapping, global_best_cost
end







using Graphs, MetaGraphsNext
includet("../utils/import_utils.jl")
includet("optimal_routing.jl")
includet("shortest-path-routing.jl")




function solve_UEPSO(instance; nb_particle=25, nb_iter=50, time_max=10, print_things=true)


    time_start = time()

    v_network = instance.v_network
    s_network = instance.s_network


    # stuff for the choice a new s s_node
    s_node_ressources = [ get_attribute_node(s_network, s_node, :cap) * sum(get_attribute_edge(s_network, get_edge(s_network, s_node, s_neighbor), :cap) for s_neighbor in neighbors(s_network, s_node)) 
                            for s_node in vertices(s_network)]

    total_ressource = sum(s_node_ressources)
    s_node_ressources = s_node_ressources / total_ressource


    position = []
    velocity = []

    personal_best = []
    personal_best_cost = []

    global_best = nothing
    global_best_cost = 9999999



    # initialization
    print_things && print("initialization... ")

    for particle in 1:nb_particle

        placement = []
        placement_cost=0
        for i in 1:nv(v_network)
            keep_on = true
            while keep_on                
                s_node = get_s_node(s_node_ressources)
                if s_node ∉ placement && get_attribute_node(s_network, s_node, :cap)>0
                    push!(placement, s_node)
                    keep_on=false
                    placement_cost+= get_attribute_node(s_network, s_node, :cost)
                end
            end
        end

        routing, routing_cost = shortest_path_routing(instance, placement)

        overall_cost = placement_cost + routing_cost

        push!(position, placement)
        push!(personal_best, position[particle])
        push!(personal_best_cost, overall_cost)

        if overall_cost < global_best_cost
            global_best = position[particle]
            global_best_cost = overall_cost
            print_things && println("We got a new best solution! value $overall_cost")
        end

        push!(velocity, ones(nv(v_network)))
    end
    print_things && println("Initialization done, best solution has cost: $global_best_cost")


    print_things && println("Starting iterations...")
    # iterations
    iter = 1
    time_total = 0
    while iter < nb_iter && time_total < time_max
        for particle in 1:nb_particle

            if personal_best_cost[particle] > 99999 # if the first isnt good, we reinitialized
                #println("We still looking...")
                placement = []
                placement_cost=0
                for i in 1:nv(v_network)
                    keep_on = true
                    while keep_on
                        s_node = get_s_node(s_node_ressources)
                        if s_node ∉ placement && get_attribute_node(s_network, s_node, :cap)>0
                            push!(placement, s_node)
                            keep_on=false
                            placement_cost+= get_attribute_node(s_network, s_node, :cost)
                        end
                    end
                end
        
                routing, routing_cost = shortest_path_routing(instance, placement)
                overall_cost = placement_cost + routing_cost
                
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
                position[particle], placement_cost = times(position[particle], velocity[particle], instance, s_node_ressources)

                routing, routing_cost = shortest_path_routing(instance, position[particle])

                overall_cost = placement_cost + routing_cost

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
        end

        iter += 1
        time_total = time() - time_start

    end

    print_things && println("UEPSO finished at iteration $nb_iter, in $(time()-time_start)s, best solution: $global_best_cost")

    if isnothing(global_best) # Need to correct this someday...
        return nothing, 99999
    end

    routing, routing_cost_shortest_path = shortest_path_routing(instance, global_best)
    final_mapping = Mapping(v_network, s_network, global_best, routing)

    return final_mapping, global_best_cost
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

    p_inertia = 0.1
    p_attraction_personal = 0.2
    p_attraction_global = 0.7

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


function times(position, velocity, instance, s_node_ressources)

    new_placement = []
    placement_cost = 0

    for i in 1:nv(instance.v_network)
        if velocity[i] == 1
            push!(new_placement, position[i])
            placement_cost += get_attribute_node(instance.s_network, position[i], :cost)
        else
            push!(new_placement, -1)
        end
    end

    for i in 1:nv(instance.v_network)
        if new_placement[i] == -1
            keep_on = true
            while keep_on
                s_node = get_s_node(s_node_ressources)
                if s_node ∉ new_placement && get_attribute_node(instance.s_network, s_node, :cap)>0
                    new_placement[i]=s_node
                    keep_on=false
                    placement_cost+= get_attribute_node(instance.s_network, s_node, :cost)
                end
            end
        end
    end

    return new_placement, placement_cost
end


function get_s_node(s_node_ressources)
    seuil = rand()

    cumul = 0.0
    for (i, val) in enumerate(s_node_ressources)
        cumul += val
        if cumul > seuil
            return i
        end
    end
    return 1
end

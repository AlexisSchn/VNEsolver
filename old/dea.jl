
using Graphs, MetaGraphsNext
includet("../utils/import_utils.jl")
includet("shortest-path-routing.jl")
includet("../utils/kahip_wrapper.jl")
includet("uepso.jl")


struct Subgraph
    graph
    real_nodes_subgraph
    neighbors_nodes_of_subgraph # nodes that are not in the subgraph, they are connected to it. But we need to put them in here for simplicity, because the edges but be done also.
    all_v_nodes
    assignment_main_to_subgraph
end



# Based on ODEA, but without overlapping decomposition.
function solve_dea(instance)

    v_network = instance.v_network
    s_network = instance.s_network

    mapping_init, cost_init = solve_UEPSO(instance)

    if cost_init > 999999
        println("Can't find an initial mapping. The road ends here...")
        return
    end

    best_placement_cost = cost_init
    best_mapping = mapping_init

    println("After initialization, we have a solution of cost $cost_init. Let's improve it...")
    # virtual graph partition
    nb_cluster = round(Int, nv(v_network)/8)
    print("Doing a partition into $nb_cluster subgraphs...")
    vn_partition = partition_vn(instance, nb_cluster)

    vn_subgraphs = []
    for (i_subgraph, v_nodes_of_subgraph) in enumerate(vn_partition)
        neighbors_of_subgraph = []
        for v_node in v_nodes_of_subgraph
            for v_neighbor in neighbors(v_network, v_node)
                if v_neighbor ∉ v_nodes_of_subgraph && v_neighbor ∉ neighbors_of_subgraph
                    push!(neighbors_of_subgraph, v_neighbor)
                end
            end 
        end
        all_v_nodes =Vector{Int64}() 
        append!(all_v_nodes, v_nodes_of_subgraph)
        append!(all_v_nodes, neighbors_of_subgraph)
        assignment_main_to_subgraph = Dict()
        for (i_v_node, v_node) in enumerate(all_v_nodes)
            assignment_main_to_subgraph[v_node] = i_v_node
        end
        println("Alll vnodes : $all_v_nodes")
        subgraph = Subgraph(my_induced_subgraph(v_network, all_v_nodes, "subgraph_$i_subgraph"), v_nodes_of_subgraph, neighbors_of_subgraph, all_v_nodes, assignment_main_to_subgraph)
        push!(vn_subgraphs, subgraph)
    end


    # let's gongue
    nb_iter = 50
    for i in 1:nb_iter
        print("iter $i... ")
        for v_subgraph in vn_subgraphs
       
            # residual graph
            residual_s_network = copy(s_network)
            for v_node in vertices(v_network)  
                if v_node ∉ v_subgraph.real_nodes_subgraph && v_node ∉ v_subgraph.real_nodes_subgraph
                    set_attribute_node(residual_s_network, v_node, :cap, 0)
                end
            end
            for v_edge in edges(v_network) 
                if src(v_edge) ∉ v_subgraph.real_nodes_subgraph && dst(v_edge) ∉ v_subgraph.real_nodes_subgraph
                    for s_edge in best_mapping.edge_routing[v_edge].edges
                        set_attribute_edge(residual_s_network, s_edge, :cap, get_attribute_edge(residual_s_network, s_edge, :cap) -1)
                    end
                end
            end
            
            # necessary placement of the neighbors nodes
            obligatory_placement = Dict()
            for v_neighbor in v_subgraph.neighbors_nodes_of_subgraph
                obligatory_placement[v_subgraph.assignment_main_to_subgraph[v_neighbor]] = best_mapping.node_placement[v_neighbor]
            end

            # Getting the submapping

            solve_submapping_UEPSO(v_subgraph.graph, residual_s_network, obligatory_placement)
            println("Well; been there, done that.")
            return


            # Integreting it to the solution, and seing if its the best?




        end

    end
end


function partition_vn(instance, nb_clusters)

    graph = instance.v_network.graph
    
    # 1 : Partitionner
    inbalance = 0.1
    println("$nb_clusters clusters to do, with inbalance $inbalance...")
    partition = partition_kahip(graph, nb_clusters, inbalance)
    clusters = [Vector{Int64}() for i in 1:nb_clusters]
    for s_node in vertices(graph)
        push!(clusters[partition[s_node]], s_node)
    end
    # 2 : Corriger
    #return clusters
    for cluster in clusters
        simple_subgraph, vmap = induced_subgraph(graph, cluster)
        if !is_connected(simple_subgraph)
            components = connected_components(simple_subgraph)
            component_sorted = sort(components, by=x->length(x), rev=true)
            for subcluster in component_sorted[2:length(component_sorted)]
                #Let's add all those nodes to a (most) connected subgraph
                #println("Look at my subcluster: $subcluster")
                nodes_original = [vmap[node] for node in subcluster]
                subgraph_neighbors = zeros(Int, nb_clusters)
                for node in nodes_original
                    for neighbor in neighbors(graph, node)
                        if neighbor ∉ cluster
                            subgraph_neighbors[partition[neighbor]] += 1
                        end
                    end
                end
                #print(subgraph_neighbors)
                most_connected_subgraph = sortperm(subgraph_neighbors, rev=true)
                cluster_to_put_nodes_in = most_connected_subgraph[1]
                append!(clusters[cluster_to_put_nodes_in], nodes_original)
                for node in nodes_original
                    partition[node] = cluster_to_put_nodes_in
                end
                #println("Well let's add $nodes_original to cluster $(clusters[most_connected_subgraph[1]])")
                filter!(e->e∉nodes_original, cluster)
            end
        end
    end
    
    # equilibrating virtual subgraphs, to get all of them pretty close...
    target_size = nv(instance.v_network) / nb_clusters
    max_size = maximum(length(cluster) for cluster in clusters)
    min_size = minimum(length(cluster) for cluster in clusters)
    iter = 1
    while max_size > 1.2 * target_size || min_size < 0.8 * target_size

        #println("Well I have some balancing to do...")
        if max_size > 1.2 * target_size # let's remove a node from the biggest subgraph
            oversized_clusters = findall(c -> length(c) > 1.2 * target_size, clusters)
            cluster_to_reduce = clusters[oversized_clusters[1]]
            print(cluster_to_reduce)

            neighboring_clusters = Dict()
            for node in cluster_to_reduce
                for neighbor in neighbors(graph, node)
                    if neighbor ∉ cluster_to_reduce
                        if partition[neighbor] ∈ keys(neighboring_clusters)
                            push!(neighboring_clusters[partition[neighbor]], node)
                        else
                            neighboring_clusters[partition[neighbor]] = [node]
                        end
                    end
                end
            end
            #println("hum sooo... $neighboring_clusters")
            smallest_cluster_around = argmin(i -> length(clusters[i]), keys(neighboring_clusters))
            #println("Well, let's add it to $smallest_cluster_around !")
            node_to_change_cluster = neighboring_clusters[smallest_cluster_around][1]
            deleteat!(cluster_to_reduce, findfirst(==(node_to_change_cluster), cluster_to_reduce)) 
            push!(clusters[smallest_cluster_around], node_to_change_cluster)
            partition[node_to_change_cluster] = smallest_cluster_around
            # put a node in the neighbor cluster with fewer nodes

        else # let's add a node to the smallest subgraph
            # do stuff..
        end

        # Update max size
        max_size = maximum(length(cluster) for cluster in clusters)
        min_size = minimum(length(cluster) for cluster in clusters)

        if iter == 1
            println("Have to do extra balancing of the subgraphs, due to KaHIP...")
        end
        iter += 1
        if iter > 30
            println("Coudn't balance the subgraphs well... Sorry.")
            break
        end
    end

    return clusters

end



# ---- heuristic for init: anything, anyone. You should be looking other things.

# ---- heuristic for iteration: it's a bit different. You need to take into account residual graph, and connecting to already existing nodes.

function solve_submapping_UEPSO(v_network, s_network, obligatory_placement)



    # stuff for the choice a new s s_node
    s_node_ressources = [ get_attribute_node(s_network, s_node, :cap) * sum(get_attribute_edge(s_network, get_edge(s_network, s_node, s_neighbor), :cap) for s_neighbor in neighbors(s_network, s_node)) 
                            for s_node in vertices(s_network)]

    total_ressource = sum(s_node_ressources)
    s_node_ressources = s_node_ressources / total_ressource



    # PSO parameters
    nb_particle=30
    nb_iter=50
    time_max = 5
    
    position = []
    velocity = []

    personal_best = []
    personal_best_cost = []

    global_best = nothing
    global_best_cost = 9999999



    # initialization
    print("initialization... ")
    for particle in 1:nb_particle

        placement = []
        placement_cost=0
        for v_node in 1:nv(v_network)
            if v_node ∈  keys(obligatory_placement)
                push!(placement, obligatory_placement[v_node])
            else
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
        end

        routing, routing_cost = shortest_path_routing(instance, placement)
        overall_cost = placement_cost + routing_cost

        push!(position, placement)
        push!(personal_best, position[particle])
        push!(personal_best_cost, overall_cost)

        if overall_cost < global_best_cost
            global_best = position[particle]
            global_best_cost = overall_cost
            println("We got a new best solution! value $overall_cost")
        end

        push!(velocity, ones(nv(v_network)))
    end
    println(" done, best solution has cost: $global_best_cost")


    println("Starting iterations...")
    # iterations
    iter = 1
    time_total = 0
    while iter < nb_iter && time_total < time_max
        for particle in 1:nb_particle

            if personal_best_cost[particle] > 99999 # if the first isnt good, we reinitialized
                #println("We still looking...")
                placement = []
                placement_cost=0
                for v_node in 1:nv(v_network)
                    if v_node ∈  keys(obligatory_placement)
                        push!(placement, obligatory_placement[v_node])
                    else
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
                    println("We got a new best solution! value $overall_cost")
                end



            else # we do a normal iteration
                velocity[particle] = plus( velocity[particle], 
                                            minus(personal_best[particle], position[particle]), 
                                            minus(global_best, position[particle]))
                position[particle], placement_cost = times_submapping(position[particle], velocity[particle], instance, s_node_ressources, obligatory_placement)

                routing, routing_cost = shortest_path_routing(instance, position[particle])

                overall_cost = placement_cost + routing_cost

                if overall_cost < personal_best_cost[particle]
                    personal_best[particle] = position[particle]
                    personal_best_cost[particle] = overall_cost
                end
                if overall_cost < global_best_cost
                    global_best_cost = overall_cost
                    global_best = position[particle]
                    println("We got a new best solution! value $global_best_cost")
                end
            end
        end

        iter += 1
        time_total = time() - time_start

    end

    #println("Final best solution: $global_best")
    println("UEPSO finished, best solution: $global_best_cost")
    #routing, routing_cost_shortest_path = shortest_path_routing(instance, global_best)
    final_mapping = Mapping(v_network, s_network, global_best, routing_cost_shortest_path)

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


function times_submapping(position, velocity, instance, s_node_ressources, obligatory_placement)

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

    for v_node in 1:nv(v_network)
        if new_placement[i] == -1

            if v_node ∈  keys(obligatory_placement)
                new_placement[i] = obligatory_placement[v_node]
            else
                keep_on = true
                while keep_on                
                    s_node = get_s_node(s_node_ressources)
                    if s_node ∉ placement && get_attribute_node(s_network, s_node, :cap)>0
                        new_placement[i]=s_node
                        keep_on=false
                        placement_cost+= get_attribute_node(s_network, s_node, :cost)
                    end
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


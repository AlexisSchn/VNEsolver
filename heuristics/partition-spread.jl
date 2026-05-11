using Random
using Graphs, MetaGraphsNext
using StatsBase

includet("../utils/import_utils.jl")
includet("../utils/partition-graph.jl")


struct Subgraph
    graph
    nodes_of_main_graph
end



function print_stuff_subgraphs(original_graph, subgraphs)

    println("For $original_graph, there is $(length(subgraphs)) subgraphs:")
    for subgraph in subgraphs
        println("       $(subgraph.graph[][:name]) with $(nv(subgraph.graph)) nodes and $(ne(subgraph.graph)) edges")
    end
    
end


struct NetworkDecomposition
    subgraphs
    v_nodes_assignment
    v_edges_master
end

function set_up_decompo(instance, node_partitionning)

    vn = instance.v_network

        
    node_assignment = Dict()
    for v_node in vertices(vn)
        node_assignment[v_node] = Dict()
    end

    # getting the subgraphs and the node assignment
    # i couldnt make the base induced_graph function work so I did adapt it
    subgraphs = []
    for (i_subgraph, v_nodes) in enumerate(node_partitionning)
        subgraph = Subgraph(my_induced_subgraph(vn, v_nodes, "subgraph_$i_subgraph"), v_nodes)
        
        for (i_node, v_node) in enumerate(v_nodes)
            node_assignment[v_node][subgraph] = i_node
        end
        push!(subgraphs, subgraph)
        #println("Look at my nice graph for the nodes $v_nodes")
        #print_graph(subgraph.graph)
    end


    # finding out the master virtual edges
    v_edge_master = [] 
    for v_edge in edges(vn)
        in_master = true
        for subgraph_src in keys(node_assignment[src(v_edge)])
            for subgraph_dst in keys(node_assignment[dst(v_edge)])
                if subgraph_src == subgraph_dst
                    in_master = false
                end
            end
        end
        if in_master
            push!(v_edge_master, v_edge)
        end
    end

    

    vn_decompo = NetworkDecomposition(subgraphs, node_assignment, v_edge_master)



    return vn_decompo
end



function complete_partial_placement_subg(partial_placement, v_network, s_network, shortest_paths_subg, capacited_nodes_subg, additional_costs_routing)


    placement = copy(partial_placement)

    #println("Partial placement : $partial_placement")
    already_placed_v_nodes = filter(v_node -> placement[v_node] != 0, vertices(v_network))


    next_v_nodes = setdiff(
        union([neighbors(v_network, v) for v in already_placed_v_nodes]...),
        already_placed_v_nodes,
    )

    possible_s_nodes = filter(s_node -> s_node ∉ placement, capacited_nodes_subg)


    while !isempty(next_v_nodes)

        # Take a node of the list
        shuffle!(next_v_nodes)
        v_node = popfirst!(next_v_nodes)

        # Get neighbors already placed
        placement_neighbors = [placement[s_neigh] for s_neigh in filter(v_neighbor -> placement[v_neighbor] != 0, neighbors(v_network, v_node))]

        # Choose some nodes
        number_s_nodes = length(possible_s_nodes)
        some_s_nodes = sample(possible_s_nodes, number_s_nodes; replace=false)


        # Rank them according distance to already placed nodes and capacity
        distances = [ sum(shortest_paths_subg.dists[s_src, s_node] for s_src in placement_neighbors) for s_node in some_s_nodes]
        distances_norm = (distances .- minimum(distances)) ./ (maximum(distances) - minimum(distances) + 1e-9)

        #capacities = [ s_node_scores[s_node] for s_node in some_s_nodes]
        #capacities_norm = (capacities .- minimum(capacities)) ./ (maximum(capacities) - minimum(capacities) + 1e-9)
        ##
        #costs = [ node_costs[s_node] + additional_costs_dual[v_node][s_node] for s_node in some_s_nodes]
        #costs_norm = (costs .- minimum(costs)) ./ (maximum(costs) - minimum(costs) + 1e-9)

        cut_routing = [ additional_costs_routing[v_node][s_node] for s_node in some_s_nodes]
        cut_norm = (cut_routing .- minimum(cut_routing)) ./ (maximum(cut_routing) - minimum(cut_routing) + 1e-9)



        final_scores = 3. .* distances_norm .+ 1. * cut_norm #.+ 0. .* costs_norm .+  0. * ( 1 .-capacities_norm)

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
        placement_cost += get_attribute_node(s_network, placement[v_node], :cost) #+ additional_costs_routing[v_node][placement[v_node]] 
    end

    return placement, placement_cost


end



function shortest_path_routing_subg(v_node_placement, v_network, s_network_dir, shortest_paths_subg, distmx_subg, capacities_edges_subg)

    capacities_edges_copy = copy(capacities_edges_subg)
    #println(s_network_dir)
    #println(capacities_edges_subg)
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
                u = shortest_paths_subg.parents[s_src, v]
                push!(nodes_of_path, u)
                v = shortest_paths_subg.parents[s_src, v]
            end
            reverse!(nodes_of_path)
            cost_of_routing_current_edge = shortest_paths_subg.dists[s_src, s_dst]
        else # A smarter thing to do here, would be to use the basic path, and if it's using a removed edge, then compute the astar path. But it does not cost so much time so.
            edges_of_path = a_star(s_network_copy_dir_copy_ofgraph, s_src, s_dst, distmx_subg)
            if edges_of_path == [] # No paths found: the graph is full!
                return Dict(), 10e9
            end
            push!(nodes_of_path, src(edges_of_path[1]))
            for edge in edges_of_path
                push!(nodes_of_path, dst(edge))
                cost_of_routing_current_edge += distmx_subg[src(edge), dst(edge)]
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
        get_routing=true 
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


function shortest_path_routing_cut_edges(overall_placement, partial_routing, 
                                                            v_network, vn_decompo, s_network_dir,
                                                            base_shortest_paths, 
                                                            dstmx, 
                                                            capacities_edges)


    # First, remove the edge capacities from the thing already done...

    complete_routing = Dict()

    capacities_edges_copy = copy(capacities_edges) 
    s_network_copy_dir_copy_ofgraph = nothing
    is_still_original_s_network = true


    for v_edge in keys(partial_routing)
        path = partial_routing[v_edge]
        # Removing the capacities. If no more capacities, removing the edge.
        for s_edge in path.edges
            capacities_edges_copy[(src(s_edge)), (dst(s_edge))] -= 1
            capacities_edges_copy[(dst(s_edge)), (src(s_edge))] -= 1

            if capacities_edges_copy[(src(s_edge)), (dst(s_edge))] == 0
                if is_still_original_s_network
                    s_network_copy_dir_copy_ofgraph = deepcopy(s_network_dir.graph)
                    is_still_original_s_network = false
                end
                rem_edge!(s_network_copy_dir_copy_ofgraph, src(s_edge), dst(s_edge))
                rem_edge!(s_network_copy_dir_copy_ofgraph, dst(s_edge), src(s_edge))
            end
        end
        complete_routing[v_edge] = path
    end 



    # then, do the new routing for the cut edgessssssss
    cost_routing_cut_edges = 0
    
    for v_edge in vn_decompo.v_edges_master

        s_src = overall_placement[src(v_edge)]
        s_dst = overall_placement[dst(v_edge)]

        nodes_of_path = Int[]
        cost_of_routing_current_edge = 0

        # Get the shortest path
        if is_still_original_s_network 
            nodes_of_path = [s_dst]
            v = s_dst
            while v != s_src
                u = base_shortest_paths.parents[s_src, v]
                push!(nodes_of_path, u)
                v = base_shortest_paths.parents[s_src, v]
            end
            reverse!(nodes_of_path)
            cost_of_routing_current_edge = base_shortest_paths.dists[s_src, s_dst]
        else # A smarter thing to do here, would be to use the basic path, and if it's using a removed edge, then compute the astar path. But it does not cost so much time so.
            edges_of_path = a_star(s_network_copy_dir_copy_ofgraph, s_src, s_dst, dstmx)
            if edges_of_path == [] # No paths found: the graph is full!
                return Dict(), 10e9
            end
            push!(nodes_of_path, src(edges_of_path[1]))
            for edge in edges_of_path
                push!(nodes_of_path, dst(edge))
                cost_of_routing_current_edge += dstmx[src(edge), dst(edge)]
            end
        end
        
        cost_routing_cut_edges += cost_of_routing_current_edge

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
        get_routing=true 
        if get_routing
            edges_of_path = Edge[]
            for i_node in 1:length(nodes_of_path)-1
                current_src = nodes_of_path[i_node]
                current_dst = nodes_of_path[i_node+1]
                push!(edges_of_path, get_edge(s_network_dir, current_src, current_dst))
            end
            path = Path(s_src, s_dst, edges_of_path, cost_of_routing_current_edge)
            complete_routing[v_edge] = path
        end
    end

    return complete_routing, cost_routing_cut_edges
end




function partition_spread(instance)



    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir




    # 1 ) Decomposing the networks

    # ------ Virtual decomposition ------ #
    nb_virtual_subgraph = floor(Int, nv(v_network.graph)/10)
    v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.2)   
    vn_decompo = set_up_decompo(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs
    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")
    
    # ------ Substrate decomposition ------ #
    nb_substrate_subgraphs = nb_virtual_subgraph
    nb_substrate_subgraphs = floor(Int, nv(s_network.graph)/20)

    clusters = partition_graph(s_network.graph, nb_substrate_subgraphs; max_umbalance = 1.3)
    sn_subgraphs = []
    for (i_subgraph, cluster) in enumerate(clusters)
        induced_subg = my_induced_subgraph(s_network, cluster, "sub_sn_$i_subgraph")
        push!(sn_subgraphs, Subgraph(induced_subg, cluster))

    end
    println("Substrate network decomposition done:")
    print_stuff_subgraphs(s_network, sn_subgraphs)


    # bruh...
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network_dir)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir, s_edge, :cost)
    end
    base_shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)

    capacities_edges = Dict{Tuple{Int, Int}, Int}()
    for s_edge in edges(s_network)
        cap = s_network[src(s_edge), dst(s_edge)][:cap]
        capacities_edges[(src(s_edge),dst(s_edge))] = cap
        capacities_edges[(dst(s_edge),src(s_edge))] = cap
    end

    sub_sn_directed = []
    shortest_paths_subgs = []
    capacities_edges_subgs = []
    dstmx_all_subg = []
    for s_subgraph in sn_subgraphs
        s_subgraph_dir = generate_dir_sn(s_subgraph.graph)
        distmx_subg = zeros(Int, nv(s_subgraph.graph), nv(s_subgraph.graph))
        cap_edges_subg = Dict{Tuple{Int, Int}, Int}()
        for s_edge in edges(s_subgraph.graph)
            cap = s_subgraph.graph[src(s_edge), dst(s_edge)][:cap]
            cap_edges_subg[(src(s_edge),dst(s_edge))] = cap
            cap_edges_subg[(dst(s_edge),src(s_edge))] = cap
        end


        for s_edge in edges(s_subgraph_dir)
            distmx_subg[src(s_edge), dst(s_edge)] = get_attribute_edge(s_subgraph_dir, s_edge, :cost)
        end
        push!(shortest_paths_subgs, floyd_warshall_shortest_paths(s_subgraph_dir, distmx_subg))
        push!(sub_sn_directed, s_subgraph_dir)
        push!(capacities_edges_subgs, cap_edges_subg)
        push!(dstmx_all_subg, distmx_subg)
    end
    
    
    
    # 2 ) you affect a random sn subg to each virtual subg
    subg_affectation = shuffle(1:nb_substrate_subgraphs)
    println("Assigment: $subg_affectation")

    # 2 ) One by one, you do the submappings
    # Don't forget to take into account the cost of previously placed neighbors!
    # And imo, you should retain that routing in the final mapping.
    temporary_placement = zeros(Int, nv(v_network))


    # initial virtua node placement, which is simple (and not correct)
    for (i_subgraph, v_subgraph) in enumerate(vn_subgraphs)
        i_s_subgraph = subg_affectation[i_subgraph]
        s_subgraph = sn_subgraphs[i_s_subgraph]
        cluster = s_subgraph.nodes_of_main_graph
        capacities_score = [ get_attribute_node(s_network, s_node, :cap) * 
                sum(get_attribute_edge(s_network, get_edge(s_network, s_node, s_neighbor), :cap) for s_neighbor in neighbors(s_network, s_node)) 
                for s_node in cluster]
        best_node = cluster[findmin(capacities_score)[2]]
        for v_node in v_subgraph.nodes_of_main_graph
            temporary_placement[v_node] = best_node
        end
    end

    # The routing will get completed, for the virtual edges in the virtual subgraphs, progressively.
    # Then, the cut edges are done at the end.
    partial_routing = Dict()

    total_cost = 0

    # Doing the submapping!
    for (i_subgraph, v_subgraph) in enumerate(vn_subgraphs)
        i_s_subgraph = subg_affectation[i_subgraph]
        s_subgraph = sn_subgraphs[i_s_subgraph]

        
        # Additional cost calculation
        additional_costs_routing = []
        for v_node in vertices(v_subgraph.graph)
            current_addition_costs = zeros(nv(s_subgraph.graph))
            original_v_node = v_subgraph.nodes_of_main_graph[v_node]
            for v_edge in vn_decompo.v_edges_master
                if src(v_edge) == original_v_node
                    placement_of_dst_node = temporary_placement[dst(v_edge)]
                    for s_node in vertices(s_subgraph.graph)
                        original_s_node = s_subgraph.nodes_of_main_graph[s_node]
                        current_addition_costs[s_node] += base_shortest_paths.dists[original_s_node, placement_of_dst_node]
                    end
                end
                if dst(v_edge) == original_v_node
                    placement_of_dst_node = temporary_placement[src(v_edge)]
                    for s_node in vertices(s_subgraph.graph)
                        original_s_node = s_subgraph.nodes_of_main_graph[s_node]
                        current_addition_costs[s_node] += base_shortest_paths.dists[original_s_node, placement_of_dst_node]
                    end
                end
            end

            push!(additional_costs_routing, current_addition_costs)
        end

        # Node placement...
        # First node placement : a bit random, kinda central...
        # Construct initial mapping
        centrality_nodes = closeness_centrality(s_subgraph.graph)
        capacited_nodes_subg = [s_node for s_node in vertices(s_subgraph.graph) if s_subgraph.graph[s_node][:cap] ≥ 1]

        s_nodes_scores = [ (centrality_nodes[s_node] + 0.5 * rand() ) for s_node in capacited_nodes_subg ]
        s_node_start = capacited_nodes_subg[argmin(s_nodes_scores)]

        partial_placement = zeros(Int32, nv(v_subgraph.graph))
        most_central_v_node = argmin(closeness_centrality(v_subgraph.graph))

        partial_placement[most_central_v_node] = s_node_start   

 


        placement_subg, placement_subg_cost = complete_partial_placement_subg(partial_placement, v_subgraph.graph, s_subgraph.graph, 
                                                                    shortest_paths_subgs[i_s_subgraph], 
                                                                    capacited_nodes_subg, 
                                                                    additional_costs_routing)

        
        # After that, do the routing, only in the SN subgraph!
        routing_subg, routing_subg_cost = shortest_path_routing_subg(placement_subg, v_subgraph.graph, 
                                                                    sub_sn_directed[i_s_subgraph], 
                                                                    shortest_paths_subgs[i_s_subgraph], 
                                                                    dstmx_all_subg[i_s_subgraph], 
                                                                    capacities_edges_subgs[i_s_subgraph])

        total_cost += placement_subg_cost + routing_subg_cost
                
        println("Placement: $placement_subg, with a cost : $(placement_subg_cost+routing_subg_cost)")

        for v_node in vertices(v_subgraph.graph)
            original_v_node = v_subgraph.nodes_of_main_graph[v_node]
            original_s_node = s_subgraph.nodes_of_main_graph[placement_subg[v_node]]
            temporary_placement[original_v_node] = original_s_node
        end
        
        for v_edge_subg in edges(v_subgraph.graph)
            original_src = v_subgraph.nodes_of_main_graph[src(v_edge_subg)]
            original_dst = v_subgraph.nodes_of_main_graph[dst(v_edge_subg)]
            original_v_edge = get_edge(v_network, original_src, original_dst)
            used_edges = []
            for s_edge in routing_subg[v_edge_subg].edges
                real_s_edge = get_edge(s_network_dir, s_subgraph.nodes_of_main_graph[src(s_edge)], s_subgraph.nodes_of_main_graph[dst(s_edge)])
                push!(used_edges, real_s_edge)
            end
            partial_routing[original_v_edge] = order_path(s_network_dir, used_edges, temporary_placement[original_src], temporary_placement[original_dst]) 
        end
        
        # And then, update the overall real routing and node placement with those.

        # Gotta write the code well, better than the old times, with all the help I can get.
        # It needs to be kinda easy to do this...


    end

    #println("Right nowwww: $temporary_placement, and routing : $partial_routing")
    # 3 ) You glue it all together with the routing in the cut edges.

    println("Cost before: $total_cost")
    overall_placement = copy(temporary_placement)
    overall_routing, cut_edge_cost = shortest_path_routing_cut_edges(overall_placement, partial_routing, 
                                                        v_network, vn_decompo, s_network_dir,
                                                        base_shortest_paths, 
                                                        distmx, 
                                                        capacities_edges)

    total_cost += cut_edge_cost

    println("Final cost: $total_cost")
    # Look at if you have a solution!
       
end

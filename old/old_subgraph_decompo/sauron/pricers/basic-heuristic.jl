# Implementation of : Virtual Network Embedding Based on the Degree and Clustering Coefficient Information
# A simple two-stage heuristic, DCC. 
# routing is done by PLNE for now, but can be done by shortest path. 
# Next: do something more progressive ?

includet("../../../utils/kahip_wrapper.jl")




function pave_heuristic(instance, vn_decompo)

    s_network = instance.s_network
    v_network = instance.v_network


    # 1 =========== central substrate nodes choices: partition the graph into, let's say, 10. Take, out of each subgraph, the most central node.

    # Substrate node ranking. The formula could surely be improved ?
    s_nodes_capacities_around = [s_network[s_node][:cap] * sum(s_network[src(s_edge), dst(s_edge)][:cap] 
                        for s_edge in get_out_edges(s_network, s_node) ∪ get_in_edges(s_network, s_node))
                            for s_node in vertices(s_network)
    ]
    s_node_centrality = closeness_centrality(s_network)
    s_nodes_scores = [s_node_centrality[s_node] * s_nodes_capacities_around[s_node] for s_node in vertices(s_network)]
    s_nodes_ordered = sortperm(s_nodes_scores; rev=true)

    # snetwork partitionning
    nb_clusters=10
    inbalance = 0.10
    partition = partition_kahip(s_network.graph, nb_clusters, inbalance)
    clusters = [Vector{Int64}() for i in 1:nb_clusters]
    for s_node in vertices(s_network)
        push!(clusters[partition[s_node]], s_node)
    end
    
    s_nodes_choiced = []
    for cluster in clusters
        for s_node in s_nodes_ordered
            if s_node in cluster
                push!(s_nodes_choiced, s_node)
                break
            end
        end
    end

    # 2 ==== CONSTRUCTING MAPPINGS
    mappings = Dict()
    for v_subgraph in vn_decompo.subgraphs
        mappings[v_subgraph] = []

        # Virtual node ranking
        v_nodes_centrality = closeness_centrality(v_subgraph.graph) # here the article does a bit different (but not really useful)
        most_central_nodes = findmax(v_nodes_centrality)[2]
        distances = desopo_pape_shortest_paths(v_subgraph.graph, most_central_nodes).dists 
        v_nodes_scores = [ v_nodes_centrality[v_node] - distances[v_node] * 10 for v_node in vertices(v_subgraph.graph)]
        v_nodes_ordered = sortperm(v_nodes_scores; rev=true)


        for s_node_central in s_nodes_choiced

            
            placement = Dict()
            v_node_placed = []
            s_node_used = []
            
            # init: we place the most central v nodes on the most central s nodes
            placement[v_nodes_ordered[1]] = s_node_central
            push!(v_node_placed, v_nodes_ordered[1])
            push!(s_node_used, s_node_central)
        
        
        
            for v_node in v_nodes_ordered[2:end]
                #println("Alright lets do $v_node")
                # find the sum of the distance between a s_node and the placements of neighbors of v_node that have been placed
                s_nodes_of_neighbors = []
                for neighbor in neighbors(v_subgraph.graph, v_node)
                    if neighbor ∈ v_node_placed
                        push!(s_nodes_of_neighbors, placement[neighbor])
                    end
                end
        
                current_distances_s_nodes = [sum(desopo_pape_shortest_paths(instance.s_network, placement_neighbor).dists[s_node] for placement_neighbor in s_nodes_of_neighbors) for s_node in vertices(s_network)]
                current_s_node_scores = [s_nodes_scores[s_node] - current_distances_s_nodes[s_node] * 10 for s_node in vertices(s_network)]
        
                current_s_nodes_ranked = sortperm(current_s_node_scores; rev=true)
                keep_on = true
                idx_s_node = 1
                while keep_on
                    s_node = current_s_nodes_ranked[idx_s_node]
                    if s_node ∉ s_node_used && s_network[s_node][:cap] >= 1 
                        s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
                        s_cap_availble = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident)
                        if s_cap_availble >= degree(v_network, v_subgraph.nodes_of_main_graph[v_node])
                            placement[v_node] = s_node
                            push!(v_node_placed, v_node)
                            push!(s_node_used, s_node)
                            keep_on = false
                        end
                    end
                    idx_s_node += 1
                end
            end
        
            #println("Placement of nodes: $placement")

            edge_routing, cost = find_routing(instance, v_subgraph.graph, placement)
        

        
            if edge_routing!==nothing

                mapping = Mapping(v_subgraph.graph, instance.s_network_dir, placement, edge_routing)

                #println("Well we find a mapping of cost $cost !")
                push!(mappings[v_subgraph], [mapping, cost])

            end



        end
        println("For subgraph $(v_subgraph.graph[][:name]), we obtained $(length(mappings[v_subgraph])) comme même ! And it was $(nv(v_subgraph.graph)) nodes")
    end

    return mappings
end











function find_routing(instance, v_network, node_placement)

    model = Model(CPLEX.Optimizer)
    set_silent(model)

    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    ### Variables
    @variable(model, x[vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);



    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);




    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
        @constraint(model, sum(x[v_node, s_node] for s_node in node_placement[v_node]) == 1)
    end

    # capacity
    for s_node in vertices(s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= s_network[s_node][:cap])
    end


    ## Edges 

    # edge capacity (undirected version !)
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( v_network[src(v_edge), dst(v_edge)][:dem] * (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
                for v_edge in edges(v_network)) 
            <= 
            instance.s_network[src(s_edge), dst(s_edge)][:cap] )
    end

    # Flow conservation
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end

    # Flow departure constraint
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end

    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        #println("Infeasible or unfinished: $status")
        return nothing, 999999
    end

    # get routing
    y_values = value.(model[:y])


    edge_routing = Dict()
    for v_edge in edges(v_network)
        used_edges = []
        for s_edge in edges(s_network_dir)
            if y_values[v_edge, s_edge] > 0.99  
                push!(used_edges, s_edge)
            end
        end
        edge_routing[v_edge] = order_path(s_network_dir, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
    end


    return edge_routing, objective_value(model)
end






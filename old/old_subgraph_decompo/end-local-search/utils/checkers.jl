



# Check if the column is valid, compared to the constraint. Returns true if it is.
function check_column_cons(column, subgraph, instance)

    s_network = instance.s_network
    s_network_dir = instance.s_network_dir
    v_network = instance.v_network

    for v_node in vertices(subgraph.graph)
        if degree(subgraph.graph, v_node) < degree(v_network, subgraph.nodes_of_main_graph[v_node])
            additional_bw_needed = degree(v_network, subgraph.nodes_of_main_graph[v_node]) - degree(subgraph.graph, v_node)
            s_node = column.node_placement[v_node]
            s_edges_incident = [get_edge(s_network_dir, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network_dir[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident;init=0.0)
            usage_edges_around = 0
            for v_edge in edges(subgraph.graph)
                path = column.edge_routing[v_edge]
                for s_edge in s_edges_incident
                    if s_edge ∈ path.edges || get_reverse_edge(s_network_dir, s_edge) ∈ path.edges
                        usage_edges_around += 1
                    end
                end
            end

            if usage_edges_around + additional_bw_needed > available_bw
                println("AHHHHHHHHHHHH")
                return false
            end
        end
    end

    return true
end




function check_column_ghost(mapping, subgraph, instance)


    # All you can do is use integer programming imo??...

    model = Model(CPLEX.Optimizer)
    set_up_checker_ghost(model, instance, subgraph, mapping)


    set_silent(model)
    optimize!(model)
    status = primal_status(model)

    if status != MOI.FEASIBLE_POINT
        println("OHHHHHHHHHHHHHHHHHHHHHHHHH")
        return false
    end

    return true

end





function set_up_checker_ghost(model, instance, subgraph, mapping)

    s_network = instance.s_network
    s_network_dir = instance.s_network_dir
    original_v_network = instance.v_network

    ghost_nodes = []
    ghost_nodes_appearances = Dict() # For each ghost node, we record how many neighbor it has in the subgraph.
    # ^ useful for additional constraints ?
    ghost_edges = []
    
    for v_node in vertices(subgraph.graph)
        for v_neighbor in neighbors(original_v_network, subgraph.nodes_of_main_graph[v_node])
            if (v_neighbor ∉ ghost_nodes) && (v_neighbor ∉ subgraph.nodes_of_main_graph)
               push!(ghost_nodes, v_neighbor) 
               ghost_nodes_appearances[v_neighbor] = 0
            end

            if v_neighbor ∉ subgraph.nodes_of_main_graph
                edge_neigh = Dict()
                edge_neigh[:src] = v_node
                edge_neigh[:dst] = v_neighbor
                push!(ghost_edges, edge_neigh)
                ghost_nodes_appearances[v_neighbor] = ghost_nodes_appearances[v_neighbor] + 1
            end
        end
    end


    #### Model

    ### Variables
    @variable(model, x[v_node in vertices(subgraph.graph), s_node in vertices(s_network)], binary=true);
    @variable(model, y[v_edge in edges(subgraph.graph), s_edge in edges(s_network_dir)], binary=true);
    @variable(model, x_ghost[ghost_nodes, vertices(s_network)], binary=true);
    @variable(model, y_ghost[ghost_edges, edges(s_network_dir)], binary=true);



    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(subgraph.graph)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
    end

    # if one to one : one virtual node per substrate node
    for s_node in vertices(s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(subgraph.graph)) <= 1)
    end



    # node capacity
    for s_node in vertices(s_network)
        @constraint(model, 
            sum( x[v_node, s_node] 
                for v_node in vertices(subgraph.graph) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end


    ## Edges 
    
    # edge capacity (undirected version)
    for s_edge in edges(s_network)
        @constraint(model, 
            sum( (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]) 
                for v_edge in edges(subgraph.graph)) 
            + sum((y_ghost[v_adjacent_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y_ghost[v_adjacent_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
                for v_adjacent_edge in ghost_edges) 
            <= 
            s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    

    # Flow conservation
    for s_node in vertices(s_network)
        for v_edge in edges(subgraph.graph)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end

    # for ghosts edges!
    for s_node in vertices(s_network)
        for v_adjacent_edge in ghost_edges
            @constraint(model, 
                x[v_adjacent_edge[:src], s_node] - x_ghost[v_adjacent_edge[:dst], s_node] 
                ==
                sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end


    ## Departure cst : Node + Edge
    for s_node in vertices(s_network)
        for v_node in vertices(subgraph.graph)
            for v_edge in get_out_edges(subgraph.graph, v_node)
                @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) >= x[v_node, s_node])
            end
        end
    end

    # for ghost edges!
    for s_node in vertices(s_network)
        for v_adjacent_edge in ghost_edges
            @constraint(model, sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[v_adjacent_edge[:src], s_node])
        end
    end


    # Outgoing edges cap: think about the original virtual graph!
    for v_node in vertices(subgraph.graph)
        for s_node in vertices(s_network)
            necessary_bw = degree(instance.v_network, subgraph.nodes_of_main_graph[v_node])
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident;init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end

    # ALSO ON GHOST NODES AHHHH
    for v_neighbor in ghost_nodes
        necessary_bw = degree(original_v_network, v_neighbor)
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, x_ghost[v_neighbor, s_node] == 0)
            end 
        end
    end


    # WELL WE ARE JUST CHECKING
    for v_node in vertices(subgraph.graph)
        @constraint(model, x[v_node, mapping.node_placement[v_node]] == 1)
    end
    for v_edge in edges(subgraph.graph)
        for s_edge in mapping.edge_routing[v_edge].edges
            @constraint(model, y[v_edge, s_edge] == 1)
        end
    end
        
end





function check_column_ghost_trois(mapping, subgraph, instance)


    # All you can do is use integer programming imo??...

    model = Model(CPLEX.Optimizer)
    set_up_checker_ghost_trois(model, instance, subgraph, mapping)


    set_silent(model)
    optimize!(model)
    status = primal_status(model)

    if status != MOI.FEASIBLE_POINT
        println("BUT NOT TROIIIIIS")
        return false
    end

    return true

end


function set_up_checker_ghost_trois(model, instance, subgraph, mapping)

    s_network = instance.s_network
    s_network_dir = instance.s_network_dir
    original_v_network = instance.v_network

    ghost_nodes = []
    ghost_nodes_appearances = Dict() # For each ghost node, we record how many neighbor it has in the subgraph.
    # ^ useful for additional constraints ?
    ghost_edges = []
    
    for v_node in vertices(subgraph.graph)
        for v_neighbor in neighbors(original_v_network, subgraph.nodes_of_main_graph[v_node])
            if (v_neighbor ∉ ghost_nodes) && (v_neighbor ∉ subgraph.nodes_of_main_graph)
               push!(ghost_nodes, v_neighbor) 
               ghost_nodes_appearances[v_neighbor] = 0
            end

            if v_neighbor ∉ subgraph.nodes_of_main_graph
                edge_neigh = Dict()
                edge_neigh[:src] = v_node
                edge_neigh[:dst] = v_neighbor
                push!(ghost_edges, edge_neigh)
                ghost_nodes_appearances[v_neighbor] = ghost_nodes_appearances[v_neighbor] + 1
            end
        end
    end


    #### Model

    ### Variables
    @variable(model, x[v_node in vertices(subgraph.graph), s_node in vertices(s_network)], binary=true);
    @variable(model, y[v_edge in edges(subgraph.graph), s_edge in edges(s_network_dir)], binary=true);
    @variable(model, x_ghost[ghost_nodes, vertices(s_network)], binary=true);
    @variable(model, y_ghost[ghost_edges, edges(s_network_dir)], binary=true);



    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(subgraph.graph)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
    end

    # if one to one : one virtual node per substrate node
    for s_node in vertices(s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(subgraph.graph)) <= 1)
    end



    # node capacity
    for s_node in vertices(s_network)
        @constraint(model, 
            sum( x[v_node, s_node] 
                for v_node in vertices(subgraph.graph) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end


    ## Edges 
    
    # edge capacity (undirected version)
    for s_edge in edges(s_network)
        @constraint(model, 
            sum( (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]) 
                for v_edge in edges(subgraph.graph)) 
            + sum((y_ghost[v_adjacent_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y_ghost[v_adjacent_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
                for v_adjacent_edge in ghost_edges) 
            <= 
            s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    

    # Flow conservation
    for s_node in vertices(s_network)
        for v_edge in edges(subgraph.graph)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end

    # for ghosts edges!
    for s_node in vertices(s_network)
        for v_adjacent_edge in ghost_edges
            @constraint(model, 
                x[v_adjacent_edge[:src], s_node] - x_ghost[v_adjacent_edge[:dst], s_node] 
                ==
                sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end


    ## Departure cst : Node + Edge
    for s_node in vertices(s_network)
        for v_node in vertices(subgraph.graph)
            for v_edge in get_out_edges(subgraph.graph, v_node)
                @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) >= x[v_node, s_node])
            end
        end
    end

    # for ghost edges!
    for s_node in vertices(s_network)
        for v_adjacent_edge in ghost_edges
            @constraint(model, sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[v_adjacent_edge[:src], s_node])
        end
    end


    # Outgoing edges cap: think about the original virtual graph!
    for v_node in vertices(subgraph.graph)
        for s_node in vertices(s_network)
            necessary_bw = degree(instance.v_network, subgraph.nodes_of_main_graph[v_node])
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident;init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end

    # ALSO ON GHOST NODES AHHHH
    for v_neighbor in ghost_nodes
        necessary_bw = degree(original_v_network, v_neighbor)
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, x_ghost[v_neighbor, s_node] == 0)
            end 
        end
    end


    # WELL WE ARE JUST CHECKING
    for v_node in vertices(subgraph.graph)
        @constraint(model, x[v_node, mapping.node_placement[v_node]] == 1)
    end
    for v_edge in edges(subgraph.graph)
        for s_edge in mapping.edge_routing[v_edge].edges
            @constraint(model, y[v_edge, s_edge] == 1)
        end
    end
        
    # I DONT WANT LONG PATHS OKAY???
    for ghost_edge in ghost_edges
        @constraint(model, sum(y_ghost[ghost_edge, s_edge] for s_edge in edges(s_network_dir)) <= 3)
    end

end







function check_column_ghost_larger(mapping, subgraph, instance)


    # All you can do is use integer programming imo??...

    model = Model(CPLEX.Optimizer)
    set_up_checker_ghost_larger(model, instance, subgraph, mapping)


    set_silent(model)
    optimize!(model)
    status = primal_status(model)

    if status != MOI.FEASIBLE_POINT
        println("BUT NOT LAAAARGER")
        return false
    end

    return true

end


function set_up_checker_ghost_larger(model, instance, subgraph, mapping)

    s_network = instance.s_network
    s_network_dir = instance.s_network_dir
    original_v_network = instance.v_network

    ghost_nodes = []
    ghost_nodes_appearances = Dict() # For each ghost node, we record how many neighbor it has in the subgraph.
    # ^ useful for additional constraints ?
    ghost_edges = []
    
    for v_node in vertices(subgraph.graph)
        for v_neighbor in neighbors(original_v_network, subgraph.nodes_of_main_graph[v_node])
            if (v_neighbor ∉ ghost_nodes) && (v_neighbor ∉ subgraph.nodes_of_main_graph)
               push!(ghost_nodes, v_neighbor) 
               ghost_nodes_appearances[v_neighbor] = 0
            end

            if v_neighbor ∉ subgraph.nodes_of_main_graph
                edge_neigh = Dict()
                edge_neigh[:src] = v_node
                edge_neigh[:dst] = v_neighbor
                push!(ghost_edges, edge_neigh)
                ghost_nodes_appearances[v_neighbor] = ghost_nodes_appearances[v_neighbor] + 1
            end
        end
    end

    ghost_ghost_edges = []
    for (i_node, v_neighbor) in enumerate(ghost_nodes)
        for (j_node, v_neighbor_2) in enumerate(ghost_nodes)
            if i_node < j_node
                if has_edge(original_v_network, v_neighbor, v_neighbor_2)
                    edge_neigh = Dict()
                    edge_neigh[:src] = v_neighbor
                    edge_neigh[:dst] = v_neighbor_2
                    push!(ghost_ghost_edges, edge_neigh)
                end
            end
        end
    end

    #### Model

    ### Variables
    @variable(model, x[v_node in vertices(subgraph.graph), s_node in vertices(s_network)], binary=true);
    @variable(model, y[v_edge in edges(subgraph.graph), s_edge in edges(s_network_dir)], binary=true);
    @variable(model, x_ghost[ghost_nodes, vertices(s_network)], binary=true);
    @variable(model, y_ghost[ghost_edges, edges(s_network_dir)], binary=true);
    @variable(model, y_ghost_ghost[ghost_ghost_edges, edges(s_network_dir)], binary=true);



    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(subgraph.graph)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
    end

    # if one to one : one virtual node per substrate node
    for s_node in vertices(s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(subgraph.graph)) <= 1)
    end



    # node capacity
    for s_node in vertices(s_network)
        @constraint(model, 
            sum( x[v_node, s_node] 
                for v_node in vertices(subgraph.graph) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end


    ## Edges 
    
    # edge capacity (undirected version)
    for s_edge in edges(s_network)
        @constraint(model, 
            sum( (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]) 
                for v_edge in edges(subgraph.graph)) 
            + sum((y_ghost[v_adjacent_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y_ghost[v_adjacent_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
                for v_adjacent_edge in ghost_edges) 
            + sum((y_ghost_ghost[v_adjacent_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y_ghost_ghost[v_adjacent_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
                for v_adjacent_edge in ghost_ghost_edges) 
            <= 
            s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    

    # Flow conservation
    for s_node in vertices(s_network)
        for v_edge in edges(subgraph.graph)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end

    # for ghosts edges!
    for s_node in vertices(s_network)
        for v_adjacent_edge in ghost_edges
            @constraint(model, 
                x[v_adjacent_edge[:src], s_node] - x_ghost[v_adjacent_edge[:dst], s_node] 
                ==
                sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end

    
    # for ghost ghost edges!
    for s_node in vertices(s_network)
        for v_adjacent_edge in ghost_ghost_edges
            @constraint(model, 
                x_ghost[v_adjacent_edge[:src], s_node] - x_ghost[v_adjacent_edge[:dst], s_node] 
                ==
                sum(y_ghost_ghost[v_adjacent_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y_ghost_ghost[v_adjacent_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end


    ## Departure cst : Node + Edge
    for s_node in vertices(s_network)
        for v_node in vertices(subgraph.graph)
            for v_edge in get_out_edges(subgraph.graph, v_node)
                @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) >= x[v_node, s_node])
            end
        end
    end

    # for ghost edges!
    for s_node in vertices(s_network)
        for v_adjacent_edge in ghost_edges
            @constraint(model, sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[v_adjacent_edge[:src], s_node])
        end
    end


    # Outgoing edges cap: think about the original virtual graph!
    for v_node in vertices(subgraph.graph)
        for s_node in vertices(s_network)
            necessary_bw = degree(instance.v_network, subgraph.nodes_of_main_graph[v_node])
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident;init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end

    # ALSO ON GHOST NODES AHHHH
    for v_neighbor in ghost_nodes
        necessary_bw = degree(original_v_network, v_neighbor)
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, x_ghost[v_neighbor, s_node] == 0)
            end 
        end
    end


    # WELL WE ARE JUST CHECKING
    for v_node in vertices(subgraph.graph)
        @constraint(model, x[v_node, mapping.node_placement[v_node]] == 1)
    end
    for v_edge in edges(subgraph.graph)
        for s_edge in mapping.edge_routing[v_edge].edges
            @constraint(model, y[v_edge, s_edge] == 1)
        end
    end
        
end

using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX, Gurobi




function basic_heuristic(instance, vn_decompo, time_max)


    model = Model(CPLEX.Optimizer)
    set_up_problem(instance, vn_decompo, model)
    set_time_limit_sec(model, time_max)
    
    for subgraph in vn_decompo.subgraphs
        for column in subgraph.columns
            add_column_ip(model, vn_decompo, instance, subgraph, column)
        end
    end

    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Infeasible or unfinished: $status")
        return -999
    end
    println("Optimal solution : $(objective_value(model))")
    return objective_value(model)
end





function set_up_problem(instance, vn_decompo, model)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)
    #set_silent(model)

    
    ### Variables
    @variable(model, x[
        v_node in vn_decompo.v_nodes_master,
        s_node in vertices(s_network)], binary = true);

    
    @variable(model, y[
        v_edge in vn_decompo.v_edges_master, 
        s_edge in edges(s_network_dir)], binary=true );
    

    lambdas = Dict()
    for subgraph in vn_decompo.subgraphs
        lambdas[subgraph] = Dict()
    end
    
    

    ### Objective
    master_placement_costs = @expression(model, sum( s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_node, s_node]
        for v_node in vn_decompo.v_nodes_master for s_node in vertices(s_network) ))

    master_routing_costs = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_edge, s_edge]
        for v_edge in vn_decompo.v_edges_master for s_edge in edges(s_network_dir) ))
    
    @objective(model, Min, master_placement_costs + master_routing_costs);

    ### Constraints

    # convexity constraints
    # Equality or inequality ? It should never be worth it... Does it ?
    @constraint(
        model, 
        mapping_selec[subgraph in vn_decompo.subgraphs],
        0 >= 1
    );

    # master virtual nodes placement
    # Equality or inequality ? It should never be worth it... Does it ?
    for v_node in vn_decompo.v_nodes_master
        @constraint(
            model,
            sum( x[v_node, s_node] for s_node in vertices(s_network)) == 1 
        )
    end



    # capacity of substrate nodes
    @constraint(
        model,
        capacity_s_node[s_node in vertices(s_network)],
        sum( v_network[v_node][:dem] * x[v_node, s_node] 
            for v_node in vn_decompo.v_nodes_master ) +
        0 
        <= s_network[s_node][:cap]
    );

    

    # capacity of substrate edges
    # undirected, so both ways !
    @constraint(
        model,
        capacity_s_edge[s_edge in edges(s_network)],
        sum( v_network[src(v_edge), dst(v_edge)][:dem] * 
                (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] +  y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))] )
            for v_edge in vn_decompo.v_edges_master)
        + 0
        <= s_network[src(s_edge), dst(s_edge)][:cap]
    );


    # flow conservation constraints
    @constraint(
        model,
        flow_conservation[v_edge in vn_decompo.v_edges_master, s_node in vertices(s_network)],
        0 == 
        sum( y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node))
        - sum( y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
    );

    # Adding x variables if nodes are in the master
    for v_edge in vn_decompo.v_edges_master
        if src(v_edge) ∈ vn_decompo.v_nodes_master
            for s_node in vertices(instance.s_network)
                set_normalized_coefficient(model[:flow_conservation][v_edge, s_node], x[src(v_edge), s_node], 1)
            end
        end
        if dst(v_edge) ∈ vn_decompo.v_nodes_master
            for s_node in vertices(instance.s_network)
                set_normalized_coefficient(model[:flow_conservation][v_edge, s_node], x[dst(v_edge), s_node], -1)
            end
        end
    end


    # Departure constraints (works only because we are in one to one !)
    @constraint(                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
        model, 
        departure[v_edge in vn_decompo.v_edges_master, s_node in vertices(s_network)],
        0 
        <=
        sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node))
    )
        
    # Adding the master x variable to the departure cst
    for v_edge in vn_decompo.v_edges_master
        if src(v_edge) ∈ vn_decompo.v_nodes_master
            for s_node in vertices(s_network)
                set_normalized_coefficient(model[:departure][v_edge, s_node], x[src(v_edge), s_node], 1)
            end
        end
    end



end



function add_column_ip(model, vn_decompo, instance, subgraph, column)

    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    lambda = @variable(model, base_name = "λ_$(subgraph.graph[][:name])_$(length(subgraph.columns))", binary=true);
    set_objective_coefficient(model, lambda, column.cost)

    # convexity
    set_normalized_coefficient(model[:mapping_selec][subgraph], lambda, 1)

    # capacities
    for s_node in vertices(s_network)

        usage = 0
        for v_node in vertices(subgraph.graph)
            if column.mapping.node_placement[v_node] == s_node
                usage += subgraph.graph[v_node][:dem]
            end
        end

        set_normalized_coefficient(model[:capacity_s_node][s_node], lambda, usage)
    end

    # Undirected, so it is a bit tricky... (transforming undirected mapping to make it better ?)
    for s_edge in edges(s_network)
        usage = 0
        for v_edge in edges(subgraph.graph)
            s_edge_one = get_edge(s_network_dir, src(s_edge), dst(s_edge))
            if s_edge_one in column.mapping.edge_routing[v_edge].edges
                usage += subgraph.graph[src(v_edge), dst(v_edge)][:dem]
            end
            s_edge_two = get_edge(s_network_dir, dst(s_edge), src(s_edge))
            if s_edge_two in column.mapping.edge_routing[v_edge].edges
                usage += subgraph.graph[src(v_edge), dst(v_edge)][:dem]
            end
        end
        set_normalized_coefficient(model[:capacity_s_edge][s_edge], lambda, usage)
    end

    
    # flow conservation 
    for v_edge in vn_decompo.v_edges_master

        if subgraph in keys(vn_decompo.v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = vn_decompo.v_nodes_assignment[src(v_edge)][subgraph]

            set_normalized_coefficient(
                model[:flow_conservation][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                1)

        end

        if subgraph in keys(vn_decompo.v_nodes_assignment[dst(v_edge)])

            v_node_in_subgraph = vn_decompo.v_nodes_assignment[dst(v_edge)][subgraph]

            set_normalized_coefficient(
                model[:flow_conservation][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                -1 )

        end


    end
    
    # departure
    for v_edge in vn_decompo.v_edges_master

        if subgraph in keys(vn_decompo.v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = vn_decompo.v_nodes_assignment[src(v_edge)][subgraph]

            set_normalized_coefficient(
                model[:departure][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                1 )

        end
    end

end





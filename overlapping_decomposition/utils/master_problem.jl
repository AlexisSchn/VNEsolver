using Graphs, MetaGraphsNext
using CPLEX, JuMP



struct MasterProblem
    instance
    model
    vn_decompo
    columns
end


struct Column
    variable
    mapping
    cost
end





function set_up_master_problem(instance, vn_decompo)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    model = Model(CPLEX.Optimizer)
    set_silent(model)

    
    ### Variables    
    @variable(model, 0. <= y[
        v_edge in vn_decompo.v_edges_master, 
        s_edge in edges(s_network_dir)] <=1. );
    

    columns = Dict()
    for subgraph in vn_decompo.subgraphs
        columns[subgraph] = []
    end
    
    
    ### Objective
    master_routing_costs = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * y[v_edge, s_edge]
        for v_edge in vn_decompo.v_edges_master for s_edge in edges(s_network_dir) ))
    
    @objective(model, Min, master_routing_costs);

    ### Constraints

    # convexity constraints
    @constraint(
        model, 
        mapping_selec[subgraph in vn_decompo.subgraphs],
        0 >= 1
    );


    # overlapping constraint
    overlapping_nodes = keys(vn_decompo.overlapping_nodes)
    @constraint(
        model,
        overlapping[s_node in vertices(s_network), v_node in overlapping_nodes, i_subgraph in 1:(length(vn_decompo.overlapping_nodes[v_node]))],
        0 <= 0
    );


    # capacity of substrate nodes
    @constraint(
        model,
        capacity_s_node[s_node in vertices(s_network)],
        0 <= s_network[s_node][:cap]
    );

    

    # capacity of substrate edges
    # undirected, so both ways !
    @constraint(
        model,
        capacity_s_edge[s_edge in edges(s_network)],
        sum( (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] +  y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))] )
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


    # Departure constraints
    @constraint(
        model, 
        departure[v_edge in vn_decompo.v_edges_master, s_node in vertices(s_network)],
        0 
        <=
        sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node))
    )
        

    add_dummy_cols(vn_decompo, model)

    return MasterProblem(instance, model, vn_decompo, columns)
end



function add_column(master_problem, instance, vn_decompo, subgraph, mapping, cost)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    model = master_problem.model

    coeff_v_nodes = [1. for i in 1:nv(subgraph.graph)]
    for v_node in vertices(subgraph.graph)
        v_node_original = subgraph.nodes_of_main_graph[v_node]
        if v_node_original in keys(vn_decompo.overlapping_nodes)
            coeff_v_nodes[v_node] = 1/length(vn_decompo.overlapping_nodes[v_node_original])
        end
    end


    lambda = @variable(model, base_name = "λ_$(subgraph.graph[][:name])_$(length(master_problem.columns[subgraph]))",lower_bound = 0., upper_bound = 1.0);
    

    real_cost = 0


    coeff_v_nodes = [1. for i in 1:nv(subgraph.graph)]
    for v_node in vertices(subgraph.graph)
        v_node_original = subgraph.nodes_of_main_graph[v_node]
        if v_node_original in keys(vn_decompo.overlapping_nodes)
            coeff_v_nodes[v_node] = 1/length(vn_decompo.overlapping_nodes[v_node_original])
        end
    end

    for v_node in vertices(subgraph.graph)
        real_cost += s_network[mapping.node_placement[v_node]][:cost] * coeff_v_nodes[v_node]
    end

    for v_edge in edges(subgraph.graph)
        real_cost += mapping.edge_routing[v_edge].cost
    end
    
    column = Column(lambda, mapping, real_cost)

    push!(master_problem.columns[subgraph], column)
    set_objective_coefficient(model, lambda, real_cost)

    # convexity
    set_normalized_coefficient(model[:mapping_selec][subgraph], lambda, 1)


    # overlapping - this is quite disgusting, but I don't see a way to do it easily...
    for v_node_overlapping in keys(vn_decompo.overlapping_nodes)
        if subgraph in vn_decompo.overlapping_nodes[v_node_overlapping]
            v_node_in_subgraph = vn_decompo.v_nodes_assignment[v_node_overlapping][subgraph]
            s_node_placement = column.mapping.node_placement[v_node_in_subgraph]

            for (i_subgraph_overlapping, subgraph_overlapping) in enumerate(vn_decompo.overlapping_nodes[v_node_overlapping])

                if subgraph_overlapping == subgraph
                    # rhs
                    set_normalized_coefficient(model[:overlapping][s_node_placement, v_node_overlapping, i_subgraph_overlapping], lambda, -1)

                    #lhs
                    if i_subgraph_overlapping == length(vn_decompo.overlapping_nodes[v_node_overlapping])
                        set_normalized_coefficient(model[:overlapping][s_node_placement, v_node_overlapping, 1], lambda, 1)
                    else
                        set_normalized_coefficient(model[:overlapping][s_node_placement, v_node_overlapping, i_subgraph_overlapping+1], lambda, 1)
                    end
                end
            end
        end
    end


    # capacities on nodes
    for s_node in vertices(s_network)
        usage = 0
        for v_node in vertices(subgraph.graph)
            if column.mapping.node_placement[v_node] == s_node
                original_v_node = subgraph.nodes_of_main_graph[v_node]
                if original_v_node in keys(vn_decompo.overlapping_nodes)
                    usage += 1/length(vn_decompo.overlapping_nodes[original_v_node])
                else
                    usage += 1
                end
            end
        end
        set_normalized_coefficient(model[:capacity_s_node][s_node], lambda, usage)
    end


    # capacities on edges (undirected case!)
    for s_edge in edges(s_network)
        usage = 0
        for v_edge in edges(subgraph.graph)
            s_edge_one = get_edge(s_network_dir, src(s_edge), dst(s_edge))
            if s_edge_one in column.mapping.edge_routing[v_edge].edges
                usage += 1
            end
            s_edge_two = get_edge(s_network_dir, dst(s_edge), src(s_edge))
            if s_edge_two in column.mapping.edge_routing[v_edge].edges
                usage += 1
            end
        end
        set_normalized_coefficient(model[:capacity_s_edge][s_edge], lambda, usage)
    end

    
    # flow conservation constraints (for vnodes that are ends of a cut edge)
    for v_edge in master_problem.vn_decompo.v_edges_master

        if subgraph in keys(master_problem.vn_decompo.v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = master_problem.vn_decompo.v_nodes_assignment[src(v_edge)][subgraph]
            if src(v_edge) in keys(vn_decompo.overlapping_nodes)
                coeff = 1/(length(vn_decompo.overlapping_nodes[src(v_edge)]))
            else
                coeff = 1
            end
            set_normalized_coefficient(
                master_problem.model[:flow_conservation][v_edge, mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                coeff)

        end

        if subgraph in keys(master_problem.vn_decompo.v_nodes_assignment[dst(v_edge)])

            v_node_in_subgraph = master_problem.vn_decompo.v_nodes_assignment[dst(v_edge)][subgraph]
            if dst(v_edge) in keys(vn_decompo.overlapping_nodes)
                coeff = 1/(length(vn_decompo.overlapping_nodes[dst(v_edge)]))
            else
                coeff = 1
            end
            set_normalized_coefficient(
                master_problem.model[:flow_conservation][v_edge, mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                -coeff )

        end


    end
    

    # departure constraints (for virtual cut edges)
    for v_edge in master_problem.vn_decompo.v_edges_master

        if subgraph in keys(master_problem.vn_decompo.v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = master_problem.vn_decompo.v_nodes_assignment[src(v_edge)][subgraph]

            if src(v_edge) in keys(vn_decompo.overlapping_nodes)
                coeff = 1/(length(vn_decompo.overlapping_nodes[src(v_edge)]))
            else
                coeff = 1
            end

            set_normalized_coefficient(
                model[:departure][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                coeff )

        end
    end

end



# DUAL COSTS!
struct DualCosts
    convexity
    overlapping
    capacity_s_node
    capacity_s_edge
    flow_conservation
    departure
end



function get_duals(instance, vn_decompo, master_problem)
    
    convexity = Dict()
    overlapping=Dict()
    capacity_s_node = Dict()
    capacity_s_edge = Dict()
    flow_conservation = Dict()
    departure = Dict()

    s_network = instance.s_network
    v_network = instance.v_network
    model = master_problem.model

    for subgraph in vn_decompo.subgraphs
        convexity[subgraph] = dual(model[:mapping_selec][subgraph])
    end

    for s_node in vertices(s_network)
        overlapping[s_node]=Dict()
        for v_node_overlapping in keys(vn_decompo.overlapping_nodes)
            overlapping[s_node][v_node_overlapping] = Dict()
            for i_subgraph in 1:length(vn_decompo.overlapping_nodes[v_node_overlapping])
                overlapping[s_node][v_node_overlapping][i_subgraph]=dual(model[:overlapping][s_node, v_node_overlapping, i_subgraph])
            end
        end
    end

    # C'est en train de se livrer à ma droite là jpp, ça parle d'amour, de mariage, de commencer relancer une relation
    # jme régale des enseignements des quarantenaires
    # - je suis grave sensible aux signes de la vie moi tu sais
    
    flow_conservation= Dict()
    #println("Flow conservation:")
    for v_edge in vn_decompo.v_edges_master
        flow_conservation[v_edge] = Dict()
        #println("   $v_edge")
        for s_node in vertices(instance.s_network)
            flow_conservation[v_edge][s_node] = dual(model[:flow_conservation][v_edge, s_node])
            #println("       $s_node: $(dual(model[:flow_conservation][v_edge, s_node]))")
        end
    end

    departure = Dict()
    #println("Departure:")
    for v_edge in vn_decompo.v_edges_master
        #println("   $v_edge")
        departure[v_edge] = Dict()
        for s_node in vertices(s_network)
            departure[v_edge][s_node] = dual(model[:departure][v_edge, s_node])
            #println("       $s_node: $(dual(model[:departure][v_edge, s_node]))")
        end
    end


    for s_node in vertices(instance.s_network)
        capacity_s_node[s_node]  = dual(master_problem.model[:capacity_s_node][s_node])
    end

    for s_edge in edges(instance.s_network)
        capacity_s_edge[s_edge]  = dual(master_problem.model[:capacity_s_edge][s_edge])
    end

    return DualCosts(convexity, overlapping, capacity_s_node, capacity_s_edge, flow_conservation, departure)
end


function get_empty_duals(instance, vn_decompo, master_problem)
    
    convexity = Dict()
    overlapping=Dict()
    capacity_s_node = Dict()
    capacity_s_edge = Dict()
    flow_conservation = Dict()
    departure = Dict()

    s_network = instance.s_network
    v_network = instance.v_network
    model = master_problem.model

    for subgraph in vn_decompo.subgraphs
        convexity[subgraph] = 0.
    end

    for s_node in vertices(s_network)
        overlapping[s_node]=Dict()
        for v_node_overlapping in keys(vn_decompo.overlapping_nodes)
            overlapping[s_node][v_node_overlapping] = Dict()
            for i_subgraph in 1:length(vn_decompo.overlapping_nodes[v_node_overlapping])
                overlapping[s_node][v_node_overlapping][i_subgraph]=0.
            end
        end
    end

    
    flow_conservation= Dict()
    #println("Flow conservation:")
    for v_edge in vn_decompo.v_edges_master
        flow_conservation[v_edge] = Dict()
        #println("   $v_edge")
        for s_node in vertices(instance.s_network)
            flow_conservation[v_edge][s_node] = 0.
            #println("       $s_node: $(dual(model[:flow_conservation][v_edge, s_node]))")
        end
    end

    departure = Dict()
    #println("Departure:")
    for v_edge in vn_decompo.v_edges_master
        #println("   $v_edge")
        departure[v_edge] = Dict()
        for s_node in vertices(s_network)
            departure[v_edge][s_node] = 0.
            #println("       $s_node: $(dual(model[:departure][v_edge, s_node]))")
        end
    end


    for s_node in vertices(instance.s_network)
        capacity_s_node[s_node]  = 0.
    end

    for s_edge in edges(instance.s_network)
        capacity_s_edge[s_edge]  = 0.
    end

    return DualCosts(convexity, overlapping, capacity_s_node, capacity_s_edge, flow_conservation, departure)
end



function average_dual_costs(instance, vn_decompo, old_dual_costs, current_dual_costs;alpha=0.8)

    convexity = Dict()
    overlapping = Dict()
    capacity_s_node = Dict()
    capacity_s_edge = Dict()
    flow_conservation = Dict()
    departure = Dict()

    s_network = instance.s_network
    v_network = instance.v_network

    for subgraph in vn_decompo.subgraphs
        convexity[subgraph] = alpha * old_dual_costs.convexity[subgraph] + (1-alpha) * current_dual_costs.convexity[subgraph]
    end


    for s_node in vertices(s_network)
        overlapping[s_node]=Dict()
        for v_node_overlapping in keys(vn_decompo.overlapping_nodes)
            overlapping[s_node][v_node_overlapping] = Dict()
            for i_subgraph in 1:length(vn_decompo.overlapping_nodes[v_node_overlapping])
                overlapping[s_node][v_node_overlapping][i_subgraph]= alpha * old_dual_costs.overlapping[s_node][v_node_overlapping][i_subgraph] + (1-alpha) * current_dual_costs.overlapping[s_node][v_node_overlapping][i_subgraph]
            end
        end
    end


    for v_edge in vn_decompo.v_edges_master
        flow_conservation[v_edge] = Dict()
        for s_node in vertices(instance.s_network)
            flow_conservation[v_edge][s_node] = alpha * old_dual_costs.flow_conservation[v_edge][s_node] + (1-alpha) * current_dual_costs.flow_conservation[v_edge][s_node]
        end
    end

    for v_edge in vn_decompo.v_edges_master
        departure[v_edge] = Dict()
        for s_node in vertices(s_network)
            departure[v_edge][s_node] = alpha * old_dual_costs.departure[v_edge][s_node]  + (1-alpha) * current_dual_costs.departure[v_edge][s_node] 
        end
    end


    for s_node in vertices(instance.s_network)
        capacity_s_node[s_node]  = alpha * old_dual_costs.capacity_s_node[s_node] + (1-alpha) * current_dual_costs.capacity_s_node[s_node]
    end

    for s_edge in edges(instance.s_network)
        capacity_s_edge[s_edge]  = alpha * old_dual_costs.capacity_s_edge[s_edge] + (1-alpha) * current_dual_costs.capacity_s_edge[s_edge]
    end


    return DualCosts(convexity, overlapping, capacity_s_node, capacity_s_edge, flow_conservation, departure)



end


### =========== INITIALIZATION of CG: dummy cols
function add_dummy_cols(vn_decompo, model)

    for subgraph in vn_decompo.subgraphs
        lambda = @variable(model, base_name = "λ_$(subgraph.graph[][:name])_dummy",lower_bound = 0., upper_bound = 1.0);
        set_objective_coefficient(model, lambda, 9999999)
    
        # constraints: only convexity!
        set_normalized_coefficient(model[:mapping_selec][subgraph], lambda, 1)
    end
end

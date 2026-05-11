
using JuMP, CPLEX

### ALL THE STUFF THAT IS NEEDED FOR THE DECOMPOSITION

### === Structs
struct NetworkDecomposition
    subgraphs
    v_nodes_assignment
    v_edges_master
end

struct Subgraph
    graph
    nodes_of_main_graph
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




############============== MASTER PROBLEM 


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


struct DualCosts
    convexity
    capacity_s_node
    capacity_s_edge
    flow_conservation
    departure
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



function add_column(master_problem, instance, subgraph, mapping, cost)

    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    model = master_problem.model

    lambda = @variable(model, base_name = "λ_$(subgraph.graph[][:name])_$(length(master_problem.columns[subgraph]))",lower_bound = 0., upper_bound = 1.0);
    column = Column(lambda, mapping, cost)

    push!(master_problem.columns[subgraph], column)
    set_objective_coefficient(model, lambda, cost)

    # convexity
    set_normalized_coefficient(model[:mapping_selec][subgraph], lambda, 1)

    # capacities on nodes
    for s_node in vertices(s_network)

        usage = 0
        for v_node in vertices(subgraph.graph)
            if column.mapping.node_placement[v_node] == s_node
                usage += 1
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

            set_normalized_coefficient(
                master_problem.model[:flow_conservation][v_edge, mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                1)

        end

        if subgraph in keys(master_problem.vn_decompo.v_nodes_assignment[dst(v_edge)])

            v_node_in_subgraph = master_problem.vn_decompo.v_nodes_assignment[dst(v_edge)][subgraph]

            set_normalized_coefficient(
                master_problem.model[:flow_conservation][v_edge, mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                -1 )

        end


    end
    
    # departure constraints (for virtual cut edges)
    for v_edge in master_problem.vn_decompo.v_edges_master

        if subgraph in keys(master_problem.vn_decompo.v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = master_problem.vn_decompo.v_nodes_assignment[src(v_edge)][subgraph]

            set_normalized_coefficient(
                model[:departure][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                1 )

        end
    end

end



function get_duals(instance, vn_decompo, master_problem)
    
    convexity = Dict()
    capacity_s_node = Dict()
    capacity_s_edge = Dict()
    flow_conservation = Dict()
    departure = Dict()

    s_network = instance.s_network
    v_network = instance.v_network
    model = master_problem.model

    convexity= Dict()
    #println("Convexity:")
    for subgraph in vn_decompo.subgraphs
        convexity[subgraph] = dual(model[:mapping_selec][subgraph])
        #println(dual(model[:mapping_selec][subgraph]))
    end

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


    return DualCosts(convexity, capacity_s_node, capacity_s_edge, flow_conservation, departure)
end



function get_empty_duals(instance, vn_decompo)

    convexity = Dict()
    capacity_s_node = Dict()
    capacity_s_edge = Dict()
    flow_conservation = Dict()
    departure = Dict()

    s_network = instance.s_network
    v_network = instance.v_network

    for subgraph in vn_decompo.subgraphs
        convexity[subgraph] = 0
    end

    for v_edge in vn_decompo.v_edges_master
        flow_conservation[v_edge] = Dict()
        for s_node in vertices(instance.s_network)
            flow_conservation[v_edge][s_node] = 0
        end
    end

    for v_edge in vn_decompo.v_edges_master
        departure[v_edge] = Dict()
        for s_node in vertices(s_network)
            departure[v_edge][s_node] = 0
        end
    end


    for s_node in vertices(instance.s_network)
        capacity_s_node[s_node]  = 0
    end

    for s_edge in edges(instance.s_network)
        capacity_s_edge[s_edge]  = 0
    end

    return DualCosts(convexity, capacity_s_node, capacity_s_edge, flow_conservation, departure)
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





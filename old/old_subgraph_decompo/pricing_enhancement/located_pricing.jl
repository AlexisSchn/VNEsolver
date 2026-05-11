using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX, Gurobi



includet("../../../utils/import_utils.jl")
#includet("../../resolution/undirected/compact_undir.jl")




# TODO:
# Repartition of costs between subgraphs
# Edges in several subgraphs
# Get the solutions nicely, of both types


### === Structs
struct NetworkDecomposition
    subgraphs
    v_nodes_assignment
    v_nodes_master
    v_edges_master
end




struct Subgraph
    graph
    nodes_of_main_graph
    columns
end



function solve_subgraph_decompo_1vn(instance, v_node_partitionning)

    v_network = instance.v_network
    s_network = instance.s_network

    println("Starting...")

    vn_decompo = set_up_decompo(instance, v_node_partitionning)

    println("Decomposition set: ")
        println("For $v_network, there is: 
            $(length(vn_decompo.subgraphs)) subgraphs, 
            $(length(vn_decompo.v_nodes_master)) nodes in no subgraph,
            $(length(vn_decompo.v_edges_master)) edges in no subgraph")


    
    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model

    print("Master problem set... ")

    nb_columns = 0
    add_dummy_cols(instance, vn_decompo, master_problem)

    print("Initial set of columns generated... ")

    pricers = Dict()
    for subgraph in vn_decompo.subgraphs
        pricers[subgraph] = set_up_pricer(instance, subgraph)
        #println("And the pricer for $subgraph is $(pricers[subgraph])")
    end
    print("Pricers set... \n")

    #------------ GENERATION DE COLONNES !
    nb_iter_max = 300
    set_silent(model)
    time_master = 0
    time_subproblems = 0
    keep_on = true
    nb_iter = 0
    best_LG = -100000
    print("\nStarting column generation: \n")

    while keep_on && nb_iter < nb_iter_max

        print("\n")
        print("Iter " * string(nb_iter))

        time1 = time()
        optimize!(model)

        status = termination_status(model)

        if status != MOI.OPTIMAL
            println("Infeasible or unfinished: $status")
            return
        end

        dual_costs = get_duals(instance, vn_decompo, master_problem)


        time2 = time()

        nb_iter += 1
        keep_on = false

        print(", CG value : " * string( floor(objective_value(master_problem.model)* 1000) / 1000) )
        total_subpb_obj = 0


        for subgraph in vn_decompo.subgraphs
            println("For subgraph $(subgraph.graph)")
            reduced_pricer(instance, vn_decompo, subgraph, dual_costs)
            column, obj_value = update_solve_pricer(instance, vn_decompo, pricers[subgraph], dual_costs)
            total_subpb_obj += obj_value
            if obj_value < -0.0001
                keep_on = true 
                add_column(master_problem, instance, subgraph, column)
                nb_columns += 1
                push!(subgraph.columns, column)
            end
            #=
            print("For subgraph $(subgraph.graph), nodes")
            for v_node in vertices(subgraph.graph)
                print(" $(column.mapping.node_placement[v_node]), ")
            end
            print("\n")
            =#
        end
    
        
    
        time3 = time()

        #Calculating LG bound
        total_dual_value = 0
        for subgraph in vn_decompo.subgraphs
            total_dual_value += dual_costs.convexity[subgraph]
        end
        for s_node in vertices(s_network)
            total_dual_value += dual_costs.capacity_s_node[s_node] * s_network[s_node][:cap]
        end
        for s_edge in edges(s_network)
            total_dual_value += dual_costs.capacity_s_edge[s_edge] * s_network[src(s_edge), dst(s_edge)][:cap]
        end
        LG_value =  total_dual_value + total_subpb_obj

        if LG_value > best_LG
            if nb_iter > 1
                best_LG = LG_value
            end
        end
        #print(", Current Lagrang.B.: " * string(floor(LG_value * 1000 ) / 1000 ))
        #print(", Best Lagrang.B.: " * string(floor(best_LG * 1000 ) / 1000 ))


        #print(", Nb Columns: " * string(nb_columns))
        time_master += time2 - time1
        time_subproblems += time3 - time2
        print("\n")
    end

        



    print("\n==================== CG finished ==================== \nReason: ")

    if !keep_on
        println("no improving columns")
    else
        println("too many columns")
    end
    println("Time in MP: " * string(time_master) * ", time in SP: " * string(time_subproblems))
    println("$nb_iter iters, final value: $(objective_value(master_problem.model))")
    println("====================================================\n")

    #println(model[:capacity_s_edge])

    ##### RECUPERATION DES SOLUTIONS
    optimize!(model)
    x_values = value.(model[:x])
    y_values = value.(model[:y])
    
    mapping_selec = Dict()
    for subgraph in vn_decompo.subgraphs
        #println("For subgraph $(subgraph.graph):\n")
        subgraph_map = Dict()
        for column in subgraph.columns
            val = value.(master_problem.lambdas[subgraph][column])
            if val > 0.001
                subgraph_map[column.mapping] = val
                #println("       $val : \n$(column.mapping)\n")
            end
        end
        mapping_selec[subgraph] = subgraph_map
    end
    
    master_node_placement = Dict()
    for v_node in vn_decompo.v_nodes_master
        node_placement = []
        for s_node in vertices(s_network)
            push!(node_placement, x_values[v_node, s_node])
        end
        master_node_placement[v_node] = node_placement
    end


    master_edge_routing = Dict()
    for v_edge in vn_decompo.v_edges_master
        edges_used = Dict()
        for s_edge in edges(instance.s_network_dir)
            if y_values[v_edge, s_edge] > 0.001
                edges_used[s_edge] = y_values[v_edge, s_edge]
            end
        end
        master_edge_routing[v_edge] = edges_used
    end
    
    mapping_subgraph = MappingDecompoFrac(instance.v_network, instance.s_network, mapping_selec, master_node_placement, master_edge_routing)
    #println("Mapping subgraph : $(mapping_subgraph)")
    

    # transform it into a classical mapping...
    #=
    node_placement = []
    for v_node in vertices(v_network)
        if v_node ∈ vn_decompo.v_nodes_master
            println("wow")
        else
            current_node_placement = zeros(length(vertices(instance.s_network)))

            for (subgraph, v_node_in_subgraph) in vn_decompo.v_nodes_assignment[v_node]
                mappings_subgraph = mapping_subgraph.subgraph_mappings[subgraph]
                for (mapping, value) in mappings_subgraph
                    current_node_placement[mapping.node_placement[v_node_in_subgraph]] += value
                end
            end

            println("For node $v_node:")
            for s_node in vertices(instance.s_network)
                if current_node_placement[s_node] > 0.001
                    println("   on node $s_node: $(current_node_placement[s_node])")
                end
            end

            push!(node_placement, current_node_placement)
        end
    end



    for v_edge in edges(v_network)
        if v_edge ∈ vn_decompo.v_edges_master
            println("Edge $(v_edge) is in master, routing")
            for (s_edge, val) in mapping_subgraph.master_routing[v_edge]
                println("   on edge $(s_edge) : $val")
            end
        else



        end
    end


    return

    #=
    edge_routing = Dict()
    for v_edge in vn_decompo.v_edges_master
        edges_used = Dict()
        for s_edge in edges(instance.s_network_dir)
            if y_values[v_edge, s_edge] > 0.001
                edges_used[s_edge] = y_values[v_edge, s_edge]
            end
        end
        master_edge_routing[v_edge] = edges_used
    end
    =#

    ####### RESOLUTION ENTIERE
    unrelax()
    optimize!(master_problem.model)
    println("Value integer : " * string(objective_value(master_problem.model)))
    =#


    return
    
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
        subgraph = Subgraph(my_induced_subgraph(vn, v_nodes, "subgraph_$i_subgraph"), v_nodes, [])
        
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
    #println("Master edges: $(v_edge_master)")

    # finding out edges accross subgraph
    v_edges_accross_subgraph = []
    for v_edge in edges(vn)
        # Get the source and destination nodes of the edge
        src_subgraphs = keys(node_assignment[src(v_edge)])  # Subgraphs of the source node
        dst_subgraphs = keys(node_assignment[dst(v_edge)])  # Subgraphs of the destination node
        
        # Find the common subgraphs between source and destination nodes
        common_subgraphs = intersect(Set(src_subgraphs), Set(dst_subgraphs))
        
        # If there are at least 2 common subgraphs, the edge is in two subgraphs
        if length(common_subgraphs) >= 2
            push!(v_edges_accross_subgraph, v_edge)
        end
    end
    # println("Edges in several subgraph: $(v_edges_accross_subgraph)))


    # finding out the master virtual nodes and nodes accross several vn
    v_nodes_in_several_subgraphs = []
    v_node_master = []
    for v_node in vertices(vn)
        if length(node_assignment[v_node]) == 0
            push!(v_node_master, v_node)
        end
        if length(node_assignment[v_node]) > 1
            push!(v_nodes_in_several_subgraphs, v_node)
        end
    end
    

    vn_decompo = NetworkDecomposition(subgraphs, node_assignment, v_node_master, v_edge_master)



    return vn_decompo
end



############============== MASTER PROBLEM 


struct MasterProblem
    instance
    model
    vn_decompo
    lambdas
end



struct Column
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
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)
    
    
    ### Variables
    @variable(model, 0. <= x[
        v_node in vn_decompo.v_nodes_master,
        s_node in vertices(s_network)] <= 1.);

    
    @variable(model, 0. <= y[
        v_edge in vn_decompo.v_edges_master, 
        s_edge in edges(s_network_dir)] <=1. );
    

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

    return MasterProblem(instance, model, vn_decompo, lambdas)
end


function add_column(master_problem, instance, subgraph, column)

    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    model = master_problem.model

    lambda = @variable(model, base_name = "λ_$(subgraph.graph[][:name])_$(length(subgraph.columns))",lower_bound = 0., upper_bound = 1.0);
    master_problem.lambdas[subgraph][column] = lambda
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
    for v_edge in master_problem.vn_decompo.v_edges_master

        if subgraph in keys(master_problem.vn_decompo.v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = master_problem.vn_decompo.v_nodes_assignment[src(v_edge)][subgraph]

            set_normalized_coefficient(
                master_problem.model[:flow_conservation][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                1)

        end

        if subgraph in keys(master_problem.vn_decompo.v_nodes_assignment[dst(v_edge)])

            v_node_in_subgraph = master_problem.vn_decompo.v_nodes_assignment[dst(v_edge)][subgraph]

            set_normalized_coefficient(
                master_problem.model[:flow_conservation][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                -1 )

        end


    end
    
    # departure
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





#############============= PRICERRRR

struct SubProblem
    model
    s_network
    v_network
    subgraph
end


function set_up_pricer(instance, subgraph)

    s_network = instance.s_network
    s_network_dir = instance.s_network_dir
    
    #### Model
    model = Model(CPLEX.Optimizer)
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)

    ### Variables
    @variable(model, x[v_node in vertices(subgraph.graph), s_node in vertices(s_network)], binary=true);
    @variable(model, y[v_edge in edges(subgraph.graph), s_edge in edges(s_network_dir)], binary=true);


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
            sum( subgraph.graph[v_node][:dem] * x[v_node, s_node] 
                for v_node in vertices(subgraph.graph) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end


    ## Edges 
    
    # edge capacity (undirected version)
    for s_edge in edges(s_network)
        @constraint(model, 
            sum( subgraph.graph[src(v_edge), dst(v_edge)][:dem] * (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]) 
                for v_edge in edges(subgraph.graph)) 
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


    ## Departure cst : Node + Edge
    for s_node in vertices(s_network)
        for v_node in vertices(subgraph.graph)
            for v_edge in get_out_edges(subgraph.graph, v_node)
                @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) >= x[v_node, s_node])
            end
        end
    end

    return SubProblem(model, instance.s_network, instance.v_network, subgraph);
end


function update_solve_pricer(instance, vn_decompo, pricer, dual_costs)

    model = pricer.model
    subgraph = pricer.subgraph
    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    ### Objective
    placement_cost = @expression(model, 
        sum( ( s_network[s_node][:cost] - dual_costs.capacity_s_node[s_node] ) * subgraph.graph[v_node][:dem] * model[:x][v_node, s_node] 
            for v_node in vertices(subgraph.graph) for s_node in vertices(s_network) ))

    routing_cost = @expression(model, sum( 
        ( s_network[src(s_edge), dst(s_edge)][:cost] - dual_costs.capacity_s_edge[s_edge] ) 
        * subgraph.graph[src(v_edge), dst(v_edge)][:dem] * (model[:y][v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + model[:y][v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))])
                for v_edge in edges(subgraph.graph) for s_edge in edges(s_network) ))


            
    # flow conservation
    flow_conservation_cost = AffExpr(0.)

    for s_node in vertices(s_network)
        for connecting_edge in vn_decompo.v_edges_master
            if subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[src(connecting_edge)][subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    -dual_costs.flow_conservation[connecting_edge][s_node] , 
                    model[:x][v_node_subgraph, s_node])
            end
            if subgraph ∈ keys(vn_decompo.v_nodes_assignment[dst(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[dst(connecting_edge)][subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    +dual_costs.flow_conservation[connecting_edge][s_node], 
                    model[:x][v_node_subgraph, s_node])
            end
        end
    end


    # departure !
    departure_costs = AffExpr(0.)
    for s_node in vertices(s_network)
        for connecting_edge in vn_decompo.v_edges_master
            if subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[src(connecting_edge)][subgraph]
                add_to_expression!(
                    departure_costs, 
                    -dual_costs.departure[connecting_edge][s_node], 
                    model[:x][v_node_subgraph,s_node])
            end
        end
    end


    @objective(model, Min, 
            -dual_costs.convexity[subgraph]
            + placement_cost + routing_cost 
            + flow_conservation_cost 
            + departure_costs);

    set_silent(model)
    optimize!(model)


    # Get the solution
    x_values = value.(model[:x])
    y_values = value.(model[:y])
    cost_of_column = 0.

    node_placement = []
    for v_node in vertices(subgraph.graph)
        for s_node in vertices(s_network)
            if x_values[v_node, s_node] > 0.99
                append!(node_placement, s_node)
                cost_of_column += subgraph.graph[v_node][:dem] * pricer.s_network[s_node][:cost]
            end
        end
    end


    edge_routing = Dict()
    for v_edge in edges(subgraph.graph)
        if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
            edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
        end
        used_edges = []
        for s_edge in edges(s_network_dir)
            if y_values[v_edge, s_edge] > 0.99
                push!(used_edges, s_edge)
                cost_of_column += subgraph.graph[src(v_edge), dst(v_edge)][:dem] * 1 * s_network_dir[src(s_edge), dst(s_edge)][:cost]
            end
        end
        edge_routing[v_edge] = order_path(s_network_dir, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
    end
    mapping = Mapping(subgraph.graph, s_network_dir, node_placement, edge_routing)
    #println(mapping)
    #println(cost_of_column)
    column = Column(mapping, cost_of_column)

    dual_value = objective_value(model)
    #println("The price of the column is : $(cost_of_column), and the obj is : $(dual_value)")
    println("   For original pb: value: $(objective_value(model)), time $(solve_time(model))")

    return column, dual_value


end




function reduced_pricer(instance, vn_decompo, subgraph, dual_costs)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir
    
    #### Model
    model = Model(CPLEX.Optimizer)
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)

    ### Variables
    @variable(model, x[v_node in vertices(subgraph.graph), s_node in vertices(s_network)], binary=true);
    @variable(model, y[v_edge in edges(subgraph.graph), s_edge in edges(s_network_dir)], binary=true);


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
            sum( subgraph.graph[v_node][:dem] * x[v_node, s_node] 
                for v_node in vertices(subgraph.graph) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end


    ## Edges 
    
    # edge capacity (undirected version)
    for s_edge in edges(s_network)
        @constraint(model, 
            sum( subgraph.graph[src(v_edge), dst(v_edge)][:dem] * (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]) 
                for v_edge in edges(subgraph.graph)) 
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


    ## Departure cst : Node + Edge
    for s_node in vertices(s_network)
        for v_node in vertices(subgraph.graph)
            for v_edge in get_out_edges(subgraph.graph, v_node)
                @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) >= x[v_node, s_node])
            end
        end
    end


    ### ========= ADDITIONAL CONSTRAINT : Constraining the node placement for nodes that are in cut edges.
    for v_node in vertices(subgraph.graph)
        to_constrain = false

        prices_s_nodes = Dict()
        for s_node in vertices(s_network)
            price = s_network[s_node][:cost] - dual_costs.capacity_s_node[s_node]

            for connecting_edge in vn_decompo.v_edges_master
                if subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                    if v_node == vn_decompo.v_nodes_assignment[src(connecting_edge)][subgraph]
                        to_constrain = true
                        price += - dual_costs.flow_conservation[connecting_edge][s_node] - dual_costs.departure[connecting_edge][s_node]
                    end
                end
                if subgraph ∈ keys(vn_decompo.v_nodes_assignment[dst(connecting_edge)])
                    if v_node == vn_decompo.v_nodes_assignment[dst(connecting_edge)][subgraph]
                        to_constrain = true
                        price += dual_costs.flow_conservation[connecting_edge][s_node]
                    end
                end
            end
            prices_s_nodes[s_node] = price
        end

        if to_constrain

            sorted_s_nodes = sort(collect(keys(prices_s_nodes)), by=x -> prices_s_nodes[x])
            @constraint(model, sum(x[v_node, s_node] for s_node in sorted_s_nodes[1:10]) == 1)

        end
    end


    ### ========== Objective : SHOULD BE DONE ELSEWHERE BUT I'M SO LAZYYYYY
    placement_cost = @expression(model, 
    sum( ( s_network[s_node][:cost] - dual_costs.capacity_s_node[s_node] ) * subgraph.graph[v_node][:dem] * model[:x][v_node, s_node] 
        for v_node in vertices(subgraph.graph) for s_node in vertices(s_network) ))

    routing_cost = @expression(model, sum( 
        ( s_network[src(s_edge), dst(s_edge)][:cost] - dual_costs.capacity_s_edge[s_edge] ) 
        * subgraph.graph[src(v_edge), dst(v_edge)][:dem] * (model[:y][v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + model[:y][v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))])
                for v_edge in edges(subgraph.graph) for s_edge in edges(s_network) ))


            
    # flow conservation
    flow_conservation_cost = AffExpr(0.)

    for s_node in vertices(s_network)
        for connecting_edge in vn_decompo.v_edges_master
            if subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[src(connecting_edge)][subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    -dual_costs.flow_conservation[connecting_edge][s_node] , 
                    model[:x][v_node_subgraph, s_node])
            end
            if subgraph ∈ keys(vn_decompo.v_nodes_assignment[dst(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[dst(connecting_edge)][subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    +dual_costs.flow_conservation[connecting_edge][s_node], 
                    model[:x][v_node_subgraph, s_node])
            end
        end
    end


    # departure !
    departure_costs = AffExpr(0.)
    for s_node in vertices(s_network)
        for connecting_edge in vn_decompo.v_edges_master
            if subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[src(connecting_edge)][subgraph]
                add_to_expression!(
                    departure_costs, 
                    -dual_costs.departure[connecting_edge][s_node], 
                    model[:x][v_node_subgraph,s_node])
            end
        end
    end


    @objective(model, Min, 
            -dual_costs.convexity[subgraph]
            + placement_cost + routing_cost 
            + flow_conservation_cost 
            + departure_costs);

    set_silent(model)
    #println("\n\n\nOKKKK SO THIS IS THE WEIRDO ONE BUDDY")

    optimize!(model)
    if is_solved_and_feasible(model)
        println("   For reduced pb  : value: $(objective_value(model)), time $(solve_time(model))")
    end
end


### ========= Column Generation result
struct MappingDecompoFrac
    v_network
    s_network
    subgraph_mappings
    master_placement
    master_routing
end


function Base.show(io::IO, mapping::MappingDecompoFrac)

    println("Mapping subgraph decomposition:")

    overall_subgraph_costs = 0
    for (subgraph, mappings) in mapping.subgraph_mappings
        current_subgraph_costs = 0

        for (mapping, val) in mappings
            current_subgraph_costs += mapping.edge_routing_cost * val + mapping.node_placement_cost * val
            if val >0.001
                println(mapping)
            end 
        end
        println("   For subgraph $(subgraph.graph), cost: $(current_subgraph_costs)")
        overall_subgraph_costs += current_subgraph_costs
    end

    println("Overall mapping of subgraphs costs : $(overall_subgraph_costs)")


    println("Master node placement:")
    for (v_node, placement) in mapping.master_placement
        println("   For vnode $(v_node)")
        for (s_node, val) in enumerate(placement)
            if val > 0.001
                println("           $(s_node)  : $(val)")
            end
        end
    end


    println("Master edge routing:")
    overall_cost_routing= 0
    for (v_edge, routing) in mapping.master_routing
        cost_for_current_edge = 0
        println("   For edge " * string(v_edge) * ":")
        for (s_edge, val) in routing
            println("           $s_edge  : $(val)")
            cost_for_current_edge += instance.s_network[src(s_edge), dst(s_edge)][:cost] * val
        end
        println("   Overall cost $(cost_for_current_edge)")
        overall_cost_routing += cost_for_current_edge
    end
    println("Overall master routing costs: $(overall_cost_routing)")


    println("Overall solution costs: $(overall_subgraph_costs + overall_cost_routing)")

end




### =========== INITIALIZATION of CG: dummy cols


function add_dummy_cols(instance, vn_decompo, master_problem)

    model = master_problem.model

    for subgraph in vn_decompo.subgraphs
        lambda = @variable(model, base_name = "λ_$(subgraph.graph[][:name])_dummy",lower_bound = 0., upper_bound = 1.0);
        set_objective_coefficient(model, lambda, 9999)
    
        # constraints: only convexity!
        set_normalized_coefficient(model[:mapping_selec][subgraph], lambda, 1)
    end
end


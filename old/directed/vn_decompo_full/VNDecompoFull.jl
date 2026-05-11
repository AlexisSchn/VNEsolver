

module NetworkDecompositionFull

using Revise

### === Includes
using Graphs, MetaGraphsNext
using JuMP, CPLEX, Gurobi

include("../../../utils/graph.jl")
include("../../../utils/import_utils.jl")



### === Structs
struct NetworkDecomposition
    subgraphs
    v_nodes_assignment
    v_nodes_master
    v_nodes_accross_subgraph
    v_edges_master
    v_edges_accross_subgraph
end




struct Subgraph
    graph
    nodes_of_main_graph
    nodes_cost_coeff
    columns
end



### === Exports
export solve_dir_vn_decompo




function solve_dir_vn_decompo(instance, v_node_partitionning)


    ### Repartition des couts entre les subpb à faire... un jour...

    println("Starting...")

    vn_decompos = set_up_decompo(instance, v_node_partitionning)
    println("Decomposition set: ")
    for v_network in instance.v_networks
        println("For $v_network, there is: 
            $(length(vn_decompos[v_network].subgraphs)) subgraphs, 
            $(length(vn_decompos[v_network].v_nodes_accross_subgraph)) nodes in several subgraphs,
            $(length(vn_decompos[v_network].v_edges_accross_subgraph)) edges accross several subgraphs,
            $(length(vn_decompos[v_network].v_nodes_master)) nodes in no subgraph,
            $(length(vn_decompos[v_network].v_edges_master)) edges in no subgraph")
        
        if length(vn_decompos[v_network].v_edges_accross_subgraph) > 0
            println("WAIT! Edges accross subgraph. This is not supported yet.")
            return 0
        end
    end




    master_problem = set_up_master_problem(instance, vn_decompos)
    print("Master problem set... ")

    nb_columns = 0
    get_initial_set_of_columns_better(instance, vn_decompos)
    get_initial_set_of_columns(instance, vn_decompos)
    for v_network in instance.v_networks
        for subgraph in vn_decompos[v_network].subgraphs
            for column in subgraph.columns
                #println("Column: $column")
                add_column(master_problem, instance, v_network, subgraph, column)
                nb_columns += 1
            end
        end
    end
    print("Initial set of columns generated... ")


    pricers = Dict()
    for v_network in instance.v_networks
        pricers[v_network] = Dict()
        for subgraph in vn_decompos[v_network].subgraphs
            pricers[v_network][subgraph] = set_up_pricer(instance, vn_decompos, v_network, subgraph)
        end
    end
    print("Pricers set... ")


    # GENERATION DE COLONNES !
    nb_iter_max = 300
    set_silent(master_problem.model)
    time_master = 0
    time_subproblems = 0
    keep_on = true
    nb_iter = 0
    best_LG = -100000
    unrelax = nothing
    print("\nStarting column generation: \n")
    while keep_on && nb_iter < nb_iter_max

        print("Iter " * string(nb_iter))
        if nb_iter > 0
            unrelax()
        end

        time1 = time()
        unrelax = relax_integrality(master_problem.model)
        optimize!(master_problem.model)
        dual_costs = get_duals(instance, vn_decompos, master_problem)
        time2 = time()

        nb_iter += 1
        keep_on = false

        print(", CG value : " * string( floor(objective_value(master_problem.model)* 1000) / 1000) )
        total_subpb_obj = 0
        for v_network in instance.v_networks
            for (i_subgraph, subgraph) in enumerate(vn_decompos[v_network].subgraphs)
                column, obj_value = update_solve_pricer(instance, pricers[v_network][subgraph], dual_costs)
                total_subpb_obj += obj_value
                if obj_value < -0.0001
                    keep_on = true 
                    add_column(master_problem, instance, v_network, subgraph, column)
                    nb_columns += 1
                    push!(subgraph.columns, column)
                end
            end
        end
        time3 = time()
        

        #Calculating LG bound
        total_dual_value = 0
        for v_network in instance.v_networks
            for subgraph in vn_decompos[v_network].subgraphs
                total_dual_value += dual_costs.convexity[v_network][subgraph]
            end
        end
        for s_node in vertices(instance.s_network)
            total_dual_value += dual_costs.capacity_s_node[s_node] * instance.s_network[s_node][:cap]
        end
        for s_edge in edges(instance.s_network)
            total_dual_value += dual_costs.capacity_s_edge[s_edge] * instance.s_network[src(s_edge), dst(s_edge)][:cap]
        end
        LG_value =  total_dual_value + total_subpb_obj
        if LG_value > best_LG
            best_LG = LG_value
        end
        print(", Lagrangian Bound: " * string(floor(best_LG * 1000 ) / 1000 ))
        

        print(", Nb Columns: " * string(nb_columns))
        time_master += time2 - time1
        time_subproblems += time3 - time2
        print("\n")
    end



    ##### RECUPERATION DES SOLUTIONS
    if nb_iter > 0
        unrelax()
    end
    unrelax = relax_integrality(master_problem.model)
    optimize!(master_problem.model)

    println("________________ \nCG finished\nFinal value: " * string(objective_value(master_problem.model)))
    println("Time in MP: " * string(time_master) * ", time in SP: " * string(time_subproblems))
    x_values = value.(master_problem.model[:x])
    y_values = value.(master_problem.model[:y])

    
    mappings_decompo = []
    
    for v_network in instance.v_networks
        mapping_selec = Dict()
        for subgraph in vn_decompos[v_network].subgraphs
            subgraph_map = Dict()
            for column in subgraph.columns
                if value.(master_problem.lambdas[v_network][subgraph][column]) > 0.001
                    subgraph_map[column.mapping] = value.(master_problem.lambdas[v_network][subgraph][column])
                end
            end
            mapping_selec[subgraph] = subgraph_map
        end
        
        
        master_node_placement = Dict()
        for v_node in vn_decompos[v_network].v_nodes_master
            node_placement = []
            for s_node in vertices(s_network)
                push!(node_placement, x_values[v_network, v_node, s_node])
            end
            master_node_placement[v_node] = node_placement
        end

        edge_routing = Dict()
        for v_edge in vn_decompos[v_network].v_edges_master
            edges_used = Dict()
            for s_edge in edges(instance.s_network)
                if y_values[v_network, v_edge, s_edge] > 0.001
                    edges_used[s_edge] = y_values[v_network, v_edge, s_edge]
                end
            end
            edge_routing[v_edge] = edges_used
        end
        

    

        
        mapping_subgraph = MappingDecompoFrac(v_network, instance.s_network, mapping_selec, master_node_placement, edge_routing)
        println("Mapping subgraph : $(mapping_subgraph)")
    end
    
    return

    # printing:
    #println(mapping[1])

    #mappings_classics = transform_into_classical_mapping(instance, mappings)
    mappings_classics = []
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            current_node_placement = zeros(length(vertices(instance.s_network)))

            for (subgraph, v_node_in_subgraph) in vn_decompos[v_network].v_nodes_assignment[v_node]
                mappings_subgraph = mappings_decompo[1].mapping_selection[subgraph]
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


    ####### RESOLUTION ENTIERE
    unrelax()
    optimize!(master_problem.model)
    #println("Value integer : " * string(objective_value(master_problem.model)))


    return mappings_decompo

end



function set_up_decompo(instance, node_partitionning)

    vn_decompos = Dict()
    
    for (i_vn, vn) in enumerate(instance.v_networks)
        
        node_assignment = Dict()
        for v_node in vertices(vn)
            node_assignment[v_node] = Dict()
        end

        # getting the subgraphs and the node assignment
        # i couldnt make the base induced_graph function work so I did adapt it
        subgraphs = []
        for (i_subgraph, v_nodes) in enumerate(node_partitionning[i_vn])
            subgraph = Subgraph(my_induced_subgraph(vn, v_nodes, "subgraph_$i_subgraph"), v_nodes, Dict(), [])
            
            for (i_node, v_node) in enumerate(v_nodes)
                node_assignment[v_node][subgraph] = i_node
            end
            push!(subgraphs, subgraph)
            println("Look at my nice graph for the nodes $v_nodes")
            print_graph(subgraph.graph)
        end


        # finding out the master virtual edges
        # for now there should not be a virtual edges on several subgraphs
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
        # Then we define the coefficient reducers...
        # Really not nice to do it here, but I can't think of doing it before blindly for now...
        # For now, equal costs among the subgraphs
        for v_node in vertices(vn)
            nb_of_subgraph_of_node = length(node_assignment[v_node])
            for (subgraph, v_node_in_subgraph) in node_assignment[v_node]
                subgraph.nodes_cost_coeff[v_node_in_subgraph] = 1/nb_of_subgraph_of_node
            end
        end

        # SIGNALER SI EDGE IN MULTI 
        vn_decompos[vn] = NetworkDecomposition(subgraphs, node_assignment, v_node_master, v_nodes_in_several_subgraphs, v_edge_master, v_edges_accross_subgraph)

    end


    return vn_decompos
end



############============== MASTER PROBLEM 


struct MasterProblem
    instance
    model
    vn_decompos
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
    splitting
end


############  ============= MASTER PROBLEM

function set_up_master_problem(instance, vn_decompos)
    model = Model(CPLEX.Optimizer)
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)
    
    
    ### Variables
    @variable(model, x[
        v_network in instance.v_networks,
        v_node in vn_decompos[v_network].v_nodes_master,
        s_node in vertices(instance.s_network)],
        binary=true);


    @variable(model, y[
        v_network in instance.v_networks, 
        v_edge in vn_decompos[v_network].v_edges_master, 
        s_edge in edges(instance.s_network)], 
        binary=true);
    
    lambdas = Dict()
    for v_network in instance.v_networks
        lambdas[v_network] = Dict()
        for subgraph in vn_decompos[v_network].subgraphs
            lambdas[v_network][subgraph] = Dict()
        end
    end
    
    

    ### Objective
    master_placement_costs = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_network, v_node, s_node]
            for v_network in instance.v_networks for v_node in vn_decompos[v_network].v_nodes_master for s_node in vertices(instance.s_network) ))

    master_routing_costs = @expression(model, sum( instance.s_network[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge]
        for v_network in instance.v_networks for v_edge in vn_decompos[v_network].v_edges_master for s_edge in edges(instance.s_network) ))
    
    @objective(model, Min, master_placement_costs + master_routing_costs);

    ### Constraints

    # convexity constraints
    @constraint(
        model, 
        mapping_selec[v_network in instance.v_networks, subgraph in vn_decompos[v_network].subgraphs],
        0 >= 1
    );

    # master virtual nodes placement
    for v_network in instance.v_networks
        for v_node in vn_decompos[v_network].v_nodes_master
            @constraint(
                model,
                sum( x[v_network, v_node, s_node] for s_node in vertices(instance.s_network)) == 1 
            )
        end
    end



    # capacity of substrate nodes
    @constraint(
        model,
        capacity_s_node[s_node in vertices(instance.s_network)],
        sum( v_network[v_node][:dem] * x[v_network, v_node, s_node] 
            for  v_network in instance.v_networks for v_node in vn_decompos[v_network].v_nodes_master ) +
        0 
        <= instance.s_network[s_node][:cap]
    );

    

    # capacity of substrate edges
    @constraint(
        model,
        capacity_s_edge[s_edge in edges(instance.s_network)],
        sum( v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge] 
            for v_network in instance.v_networks for v_edge in vn_decompos[v_network].v_edges_master)
        + 0
        <= instance.s_network[src(s_edge), dst(s_edge)][:cap]
    );


    # flow conservation constraints
    @constraint(
        model,
        flow_conservation[v_network in instance.v_networks, v_edge in vn_decompos[v_network].v_edges_master, s_node in vertices(instance.s_network)],
        0 == 
        sum( y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node))
        - sum( y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node))
    );


    # Departure constraints (works only because we are in one to one !!!!!)
    @constraint(
        model, 
        departure[v_network in instance.v_networks, v_edge in vn_decompos[v_network].v_edges_master, s_node in vertices(instance.s_network)],
        0 
        <=
        sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node))
    )

    # Adding the master x variable to the two constraints
    for v_network in instance.v_networks
        for v_edge in vn_decompos[v_network].v_edges_master
            if src(v_edge) ∈ vn_decompos[v_network].v_nodes_master
                for s_node in vertices(instance.s_network)
                    set_normalized_coefficient(model[:flow_conservation][v_network, v_edge, s_node], x[v_network, src(v_edge), s_node], 1)
                    set_normalized_coefficient(model[:departure][v_network, v_edge, s_node], x[v_network, src(v_edge), s_node], 1)
                end
            end
            if dst(v_edge) ∈ vn_decompos[v_network].v_nodes_master
                for s_node in vertices(instance.s_network)
                    set_normalized_coefficient(model[:flow_conservation][v_network, v_edge, s_node], x[v_network, dst(v_edge), s_node], -1)
                end
            end
        end
    end

    # Variable splitting variable and constraint
    @variable(model,
        0 <= x_splitted[v_network in instance.v_networks,
        v_node in vn_decompos[v_network].v_nodes_accross_subgraph,
        s_node in vertices(instance.s_network)]
        <= 1)
    

    @constraint(
        model, 
        splitting[
            v_network in instance.v_networks, 
            v_node in vn_decompos[v_network].v_nodes_accross_subgraph, 
            s_node in vertices(instance.s_network), 
            subgraph in keys(vn_decompos[v_network].v_nodes_assignment[v_node])],
        0
        ==
        x_splitted[v_network, v_node, s_node]
    )

    #=
    @constraint(
        model,
        [v_network in instance.v_networks, 
        v_node in vn_decompos[v_network].v_nodes_accross_subgraph, 
        s_node in vertices(instance.s_network)],
        sum(x_splitted[v_network, v_node, s_node] for s_node in vertices(instance.s_network)) == 1
    )
    =#


    ## Additional constraints : Undirected edges, thus the value both way is the same
    undirected = false
    if undirected
        for v_network in instance.v_networks
            for v_edge in vn_decompos[v_network].v_edges_master
                for s_edge in edges(instance.s_network)
                    @constraint(model, y[v_network, v_edge, s_edge] == y[v_network, get_edge(v_network, dst(v_edge), src(v_edge)), 
                                                get_edge(instance.s_network, dst(s_edge), src(s_edge))])
                end
            end
        end
    end

    return MasterProblem(instance, model, vn_decompos, lambdas)
end


function add_column(master_problem, instance, v_network, subgraph, column)

    lambda = @variable(master_problem.model, binary=true);
    master_problem.lambdas[v_network][subgraph][column] = lambda
    set_objective_coefficient(master_problem.model, lambda, column.cost)

    # convexity
    set_normalized_coefficient(master_problem.model[:mapping_selec][v_network, subgraph], lambda, 1)

    
    # capacities
    for s_node in vertices(instance.s_network)

        usage = 0
        for v_node in vertices(subgraph.graph)
            if column.mapping.node_placement[v_node] == s_node
                usage += subgraph.graph[v_node][:dem] * subgraph.nodes_cost_coeff[v_node]
            end
        end

        set_normalized_coefficient(master_problem.model[:capacity_s_node][s_node], lambda, usage)
    end

    for s_edge in edges(instance.s_network)
        usage = 0
        for v_edge in edges(subgraph.graph)
            if s_edge in column.mapping.edge_routing[v_edge].edges
                usage += subgraph.graph[src(v_edge), dst(v_edge)][:dem] * 1 # no splitting of edges among several subgraphs... for now
            end
        end
        set_normalized_coefficient(master_problem.model[:capacity_s_edge][s_edge], lambda, usage)
    end

    
    # flow conservation 
    for v_edge in master_problem.vn_decompos[v_network].v_edges_master

        if subgraph in keys(master_problem.vn_decompos[v_network].v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = master_problem.vn_decompos[v_network].v_nodes_assignment[src(v_edge)][subgraph]

            set_normalized_coefficient(
                master_problem.model[:flow_conservation][v_network, v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                1*subgraph.nodes_cost_coeff[v_node_in_subgraph] )

        end

        if subgraph in keys(master_problem.vn_decompos[v_network].v_nodes_assignment[dst(v_edge)])

            v_node_in_subgraph = master_problem.vn_decompos[v_network].v_nodes_assignment[dst(v_edge)][subgraph]

            set_normalized_coefficient(
                master_problem.model[:flow_conservation][v_network, v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                -1*subgraph.nodes_cost_coeff[v_node_in_subgraph] )

        end


    end
    
    # departure
    for v_edge in master_problem.vn_decompos[v_network].v_edges_master

        if subgraph in keys(master_problem.vn_decompos[v_network].v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = master_problem.vn_decompos[v_network].v_nodes_assignment[src(v_edge)][subgraph]

            set_normalized_coefficient(
                master_problem.model[:departure][v_network, v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                1*subgraph.nodes_cost_coeff[v_node_in_subgraph] )

        end
    end


    # Variable splitting: only if the v_node is shared among vn
    # then we add it to the corresponding constraint in the master.
    for v_node in vertices(subgraph.graph)
        v_node_in_original_graph = subgraph.nodes_of_main_graph[v_node] 
        if v_node_in_original_graph in master_problem.vn_decompos[v_network].v_nodes_accross_subgraph
            set_normalized_coefficient(
                master_problem.model[:splitting][v_network, v_node_in_original_graph, column.mapping.node_placement[v_node], subgraph], 
                lambda, 
                1)

        end
    end
end



function get_duals(instance, vn_decompos, master_problem)
    
    convexity = Dict()
    flow_conservation = Dict()
    departure = Dict()
    splitting = Dict()

    for v_network in instance.v_networks

        convexity[v_network] = Dict()
        for subgraph in vn_decompos[v_network].subgraphs
            convexity[v_network][subgraph] = dual(master_problem.model[:mapping_selec][v_network, subgraph])
        end

        flow_conservation[v_network] = Dict()
        for v_edge in vn_decompos[v_network].v_edges_master
            flow_conservation[v_network][v_edge] = Dict()
            for s_node in vertices(instance.s_network)
                flow_conservation[v_network][v_edge][s_node] = dual(master_problem.model[:flow_conservation][v_network, v_edge, s_node])
            end
        end

        departure[v_network] = Dict()
        for v_edge in vn_decompos[v_network].v_edges_master
            departure[v_network][v_edge] = Dict()
            for s_node in vertices(instance.s_network)
                departure[v_network][v_edge][s_node] = dual(master_problem.model[:departure][v_network, v_edge, s_node])
            end
        end

        splitting[v_network] = Dict()
        for v_node in vn_decompos[v_network].v_nodes_accross_subgraph
            splitting[v_network][v_node] = Dict()
            for s_node in vertices(instance.s_network)
                splitting[v_network][v_node][s_node] = Dict()
                for subgraph in keys(vn_decompos[v_network].v_nodes_assignment[v_node])
                    splitting[v_network][v_node][s_node][subgraph] = dual(master_problem.model[:splitting][v_network, v_node, s_node, subgraph])
                end
            end
        end

    end
    
    capacity_s_node = Dict()

    for s_node in vertices(instance.s_network)
        capacity_s_node[s_node]  = dual(master_problem.model[:capacity_s_node][s_node])
    end

    capacity_s_edge = Dict()
    for s_edge in edges(instance.s_network)
        capacity_s_edge[s_edge]  = dual(master_problem.model[:capacity_s_edge][s_edge])
    end


    return DualCosts(convexity, capacity_s_node, capacity_s_edge, flow_conservation, departure, splitting)
end





#############============= PRICERRRR

struct SubProblem
    model
    vn_decompos
    s_network
    v_network
    subgraph
end


function set_up_pricer(instance, vn_decompos, v_network_used, subgraph)

    s_network = instance.s_network

    #### Model
    model = Model(CPLEX.Optimizer)
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)

    ### Variables
    @variable(model, x[v_node in vertices(subgraph.graph), s_node in vertices(s_network)], binary=true);
    @variable(model, y[v_edge in edges(subgraph.graph), s_edge in edges(s_network)], binary=true);


    ### Constraints


    one_to_one = true
    departure_cst = true

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(subgraph.graph)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
    end

    # if one to one : one virtual node per substrate node
    if one_to_one
        for s_node in vertices(s_network)
            @constraint(model, sum(x[v_node, s_node] for v_node in vertices(subgraph.graph)) <= 1)
        end
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
    
    # edge capacity
    for s_edge in edges(s_network)
        @constraint(model, 
            sum( subgraph.graph[src(v_edge), dst(v_edge)][:dem] * y[v_edge, s_edge] 
                for v_edge in edges(subgraph.graph)) 
            <= 
            s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(s_network)
        for v_edge in edges(subgraph.graph)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                <=
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network, s_node))
            )
        end
    end


    ## Additional constraints : Node + Edge
    if one_to_one && departure_cst
        for s_node in vertices(s_network)
            for v_node in vertices(subgraph.graph)
                for v_edge in get_out_edges(subgraph.graph, v_node)
                    @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network, s_node)) >= x[v_node, s_node])
                end
                for v_edge in get_in_edges(subgraph.graph, v_node)
                    @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network, s_node)) >= x[v_node, s_node])
                end
            end
        end
    end

    ## Additional constraints : Undirected edges, thus the value both way is the same
    undirected = true
    if undirected
        for v_edge in edges(subgraph.graph)
            for s_edge in edges(s_network)
                @constraint(model, y[v_edge, s_edge] == y[get_edge(subgraph.graph, dst(v_edge), src(v_edge)), 
                                            get_edge(s_network, dst(s_edge), src(s_edge))])
            end
        end
    end



    return SubProblem(model, vn_decompos, instance.s_network, v_network_used, subgraph);
end


function update_solve_pricer(instance, pricer, dual_costs)

    model = pricer.model
    subgraph = pricer.subgraph
    vn_decompos = pricer.vn_decompos

    ### Objective
    placement_cost = @expression(model, 
        sum( ( pricer.s_network[s_node][:cost] - dual_costs.capacity_s_node[s_node] ) * subgraph.graph[v_node][:dem] * subgraph.nodes_cost_coeff[v_node] * model[:x][v_node, s_node] 
            for v_node in vertices(subgraph.graph) for s_node in vertices(pricer.s_network) ))

    routing_cost = @expression(model, sum( ( 
                        pricer.s_network[src(s_edge), dst(s_edge)][:cost] 
                        - dual_costs.capacity_s_edge[s_edge] ) * subgraph.graph[src(v_edge), dst(v_edge)][:dem] * model[:y][v_edge, s_edge] 
        for v_edge in edges(subgraph.graph) for s_edge in edges(pricer.s_network) ))


            
    # flow conservation
    flow_conservation_cost = AffExpr(0.)

    for s_node in vertices(pricer.s_network)
        for connecting_edge in vn_decompos[pricer.v_network].v_edges_master
            if subgraph ∈ keys(vn_decompos[pricer.v_network].v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompos[pricer.v_network].v_nodes_assignment[src(connecting_edge)][subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    -dual_costs.flow_conservation[pricer.v_network][connecting_edge][s_node] * subgraph.nodes_cost_coeff[v_node_subgraph], 
                    model[:x][vn_decompos[pricer.v_network].v_nodes_assignment[src(connecting_edge)][subgraph], s_node])
            end
            if subgraph ∈ keys(vn_decompos[pricer.v_network].v_nodes_assignment[dst(connecting_edge)])
                v_node_subgraph = vn_decompos[pricer.v_network].v_nodes_assignment[dst(connecting_edge)][subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    +dual_costs.flow_conservation[pricer.v_network][connecting_edge][s_node] * subgraph.nodes_cost_coeff[v_node_subgraph], 
                    model[:x][vn_decompos[pricer.v_network].v_nodes_assignment[dst(connecting_edge)][subgraph], s_node])
            end
        end
    end

    # splitting
    splitting_costs = AffExpr(0.)
    for v_node in vn_decompos[pricer.v_network].v_nodes_accross_subgraph
        if subgraph in keys(vn_decompos[pricer.v_network].v_nodes_assignment[v_node])
            v_node_subgraph = vn_decompos[pricer.v_network].v_nodes_assignment[v_node][subgraph]
            for s_node in vertices(pricer.s_network)
                add_to_expression!(
                    splitting_costs,
                    -dual_costs.splitting[pricer.v_network][v_node][s_node][subgraph],
                    model[:x][v_node_subgraph, s_node]
                )
            end
        end
    end

    # departure !
    departure_costs = AffExpr(0.)
    for s_node in vertices(pricer.s_network)
        for connecting_edge in vn_decompos[pricer.v_network].v_edges_master
            if subgraph ∈ keys(vn_decompos[pricer.v_network].v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompos[pricer.v_network].v_nodes_assignment[src(connecting_edge)][subgraph]
                add_to_expression!(
                    departure_costs, 
                    -dual_costs.departure[pricer.v_network][connecting_edge][s_node] * subgraph.nodes_cost_coeff[v_node_subgraph], 
                    model[:x][v_node_subgraph,s_node])
            end
        end
    end


    @objective(model, Min, 
            -dual_costs.convexity[pricer.v_network][subgraph]
            + placement_cost + routing_cost 
            + flow_conservation_cost 
            + departure_costs
            + splitting_costs);


    set_silent(model)
    optimize!(model)


    # Get the solution
    x_values = value.(pricer.model[:x])
    y_values = value.(pricer.model[:y])
    cost_of_column = 0.

    node_placement = []
    for v_node in vertices(subgraph.graph)
        for s_node in vertices(pricer.s_network)
            if x_values[v_node, s_node] > 0.99
                append!(node_placement, s_node)
                cost_of_column += subgraph.graph[v_node][:dem] * subgraph.nodes_cost_coeff[v_node] * pricer.s_network[s_node][:cost]
            end
        end
    end

    edge_routing = Dict()
    for v_edge in edges(subgraph.graph)
        if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
            edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
        end
        used_edges = []
        for s_edge in edges(instance.s_network)
            if y_values[v_edge, s_edge] > 0.99
                push!(used_edges, s_edge)
                cost_of_column += subgraph.graph[src(v_edge), dst(v_edge)][:dem] * 1 * pricer.s_network[src(s_edge), dst(s_edge)][:cost]
            end
        end
        edge_routing[v_edge] = order_path(pricer.s_network, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
    end
    mapping = Mapping(subgraph.graph, instance.s_network, node_placement, edge_routing)
    

    column = Column(mapping, cost_of_column)

    dual_value = objective_value(model)
        
    return column, dual_value


end





### =========== INITIALIZATION of CG: creating the first set of columns
# This is pretty archaic, surely there is some better heuristic way to do this. But at least it works in almost every cases (i.e. if the edge capacity is large enough)

# for unit capacity only (need small adaptation, remove the used cap)
# A bit old, be careful

### Modifies vn_decompos. Shouldnt it be with ! in the function ?
function get_initial_set_of_columns(instance, vn_decompos)
    placement_so_far = Dict()
    used_nodes = []
    forced_placement = Dict()

    for v_network in instance.v_networks
        placement_so_far[v_network] = Dict()
        for subgraph in vn_decompos[v_network].subgraphs
            current_instance = InstanceVNE([subgraph.graph], instance.s_network)
            forced_placement[subgraph.graph] = Dict()
            for v_node in keys(placement_so_far[v_network])
                if subgraph ∈ keys(vn_decompos[v_network].v_nodes_assignment[v_node])
                    forced_placement[subgraph.graph][vn_decompos[v_network].v_nodes_assignment[v_node][subgraph]] = placement_so_far[v_network][v_node]
                end
            end
            cost_coeff = Dict()
            cost_coeff[subgraph.graph] = subgraph.nodes_cost_coeff
            column = solve_integer_with_forbidden_nodes_and_force_placement(current_instance, cost_coeff, used_nodes, forced_placement)[1]
            for (v_node, s_node) in enumerate(column.mapping.node_placement)
                if s_node ∉ used_nodes
                    push!(used_nodes, s_node)
                    placement_so_far[v_network][subgraph.nodes_of_main_graph[v_node]] = s_node            
                end
            end
            push!(subgraph.columns, column);
            
            print(column.mapping)
        end
    end
end


# Todo: use the base function...
function solve_integer_with_forbidden_nodes_and_force_placement(instance, nodes_cost_coeff, forbidden_nodes, forced_placement)

    
    #### Model
    model = Model(CPLEX.Optimizer)
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)

    ### Variables
    @variable(model, x[v_network in instance.v_networks, vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[v_network in instance.v_networks, edges(v_network), edges(instance.s_network)], binary=true);

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_network, v_node, s_node] * nodes_cost_coeff[v_network][v_node]
        for v_network in instance.v_networks for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( instance.s_network[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge]
        for v_network in instance.v_networks for v_edge in edges(v_network) for s_edge in edges(instance.s_network) ))
    @objective(model, Min, placement_cost + routing_cost);


    ### Constraints     

    one_to_one = true
    departure_cst = true


    ### ======== Additional constraints
    
    ### Forced placement
    for v_network in instance.v_networks
        for v_node in keys(forced_placement[v_network])
            @constraint(model, x[v_network, v_node, forced_placement[v_network][v_node]] == 1)    
        end
    end

    ### forbidden nodes (not for already placed nodes)
    for s_node in forbidden_nodes
        for v_network in instance.v_networks
            for v_node in vertices(v_network)
                if v_node ∉ keys(forced_placement[v_network])
                    @constraint(model, x[v_network, v_node, s_node] == 0)
                end
            end
        end
    end


    ### =========== Nodes constraints

    # one substrate node per virtual node
    for v_network in instance.v_networks
        for v_node in vertices(v_network)
            @constraint(model, sum(x[v_network, v_node, s_node] for s_node in vertices(instance.s_network)) == 1)
        end
    end

    # if one to one : one virtual node per substrate node
    if one_to_one
        for s_node in vertices(instance.s_network)
            for v_network in instance.v_networks
                @constraint(model, sum(x[v_network, v_node, s_node] for v_node in vertices(v_network)) <= 1)
            end
        end
    end

    # node capacity
    for s_node in vertices(instance.s_network)
        @constraint(model, 
            sum( v_network[v_node][:dem] * x[v_network, v_node, s_node] 
                for v_network in instance.v_networks for v_node in vertices(v_network) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end


    ### ========== Edges constraints 
    
    # edge capacity
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge] 
                for v_network in instance.v_networks for v_edge in edges(v_network)) 
            <= 
            instance.s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(instance.s_network)
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                @constraint(model, 
                    x[v_network, src(v_edge), s_node] - x[v_network, dst(v_edge), s_node] 
                    <=
                    sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) - 
                        sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node))
                )
            end
        end
    end


    ### ========= Additional constraints
    
    ## Departure constraint
    if one_to_one
        if departure_cst
            for s_node in vertices(instance.s_network)
                for v_network in instance.v_networks
                    for v_node in vertices(v_network)
                        for v_edge in get_out_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                        end
                        for v_edge in get_in_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                        end
                    end
                end
            end
        end
    end

    ## Symmetric edges routings (because undirected !)
    undirected = true
    if undirected
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                for s_edge in edges(instance.s_network)
                    @constraint(model, y[v_network, v_edge, s_edge] == y[v_network, get_edge(v_network, dst(v_edge), src(v_edge)), 
                                                get_edge(instance.s_network, dst(s_edge), src(s_edge))])
                end
            end
        end
    end

    # Solving
    set_silent(model)
    optimize!(model)

    # Get the solution
    x_values = value.(model[:x])
    y_values = value.(model[:y])
    columns = []
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            for s_node in vertices(instance.s_network)
                if x_values[v_network, v_node, s_node] > 0.99
                    append!(node_placement, s_node)
                end
            end
        end

        edge_routing = Dict()
        for v_edge in edges(v_network)
            if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
                edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
            end
            used_edges = []
            for s_edge in edges(instance.s_network)
                if y_values[v_network, v_edge, s_edge] > 0.99
                    push!(used_edges, s_edge)
                end
            end
            edge_routing[v_edge] = order_path(instance.s_network, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
        end
        m = Mapping(v_network, instance.s_network, node_placement, edge_routing)

        column = Column(m, objective_value(model))
        push!(columns, column)
    end

    return columns
end



### A better initialization: trying to put each vn a bit everywhere
# does not work for overlapping decompo for now
function get_initial_set_of_columns_better(instance, vn_decompos)

    nb_column_per_subgraph = 10

    nb_s_nodes = length(vertices(instance.s_network))

    #substrate_nodes = shuffle(L) # for random
    substrate_nodes = 1:nb_s_nodes

    # Compute the approximate size of each partition
    base_size = div(nb_s_nodes, nb_column_per_subgraph)
    extra = nb_s_nodes % nb_column_per_subgraph

    # Partition the shuffled list into roughly equal parts
    groups = []
    start_idx = 1
    for i in 1:nb_column_per_subgraph
        # Distribute the remainder across the first few groups
        group_size = base_size + (i <= extra ? 1 : 0)
        end_idx = start_idx + group_size - 1
        push!(groups, substrate_nodes[start_idx:end_idx])
        start_idx = end_idx + 1
    end

    #println("groups: $groups")
    # get the most central node
    # assign it somewhere on the selected nodes.
    # => nodes at random ? It would be better through some gentle clustering ?
    # No need to be just one node, it should be much better if it has some choice. I think random would be good to begin with. 
    # solve the plne, add the column.

    for v_network in instance.v_networks
        for subgraph in vn_decompos[v_network].subgraphs
            current_instance = InstanceVNE([subgraph.graph], instance.s_network)
            cost_coeff = Dict()
            cost_coeff[subgraph.graph] = subgraph.nodes_cost_coeff
            println("Let's have fun !")
            for group in groups
                # the virtual node should be not too connected to make it easy ?
                print(group)
                placement_restriction = Dict()
                placement_restriction[1] = group
                # create the plne
                column = solve_integer_with_placement_restriction(current_instance, cost_coeff, placement_restriction)
                # add the column
                
                if column !== nothing
                    push!(subgraph.columns, column);
                end
            end
        end
    end

end

# Todo: use the base function...
function solve_integer_with_placement_restriction(instance, nodes_cost_coeff, placement_restriction)

    
    #### Model
    model = Model(CPLEX.Optimizer)
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)

    ### Variables
    @variable(model, x[v_network in instance.v_networks, vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[v_network in instance.v_networks, edges(v_network), edges(instance.s_network)], binary=true);

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_network, v_node, s_node] * nodes_cost_coeff[v_network][v_node]
        for v_network in instance.v_networks for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( instance.s_network[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge]
        for v_network in instance.v_networks for v_edge in edges(v_network) for s_edge in edges(instance.s_network) ))
    @objective(model, Min, placement_cost + routing_cost);


    ### Constraints     

    one_to_one = true
    departure_cst = true


    ### ======== Additional constraints
    ### placement_restriction
    for v_network in instance.v_networks
        for v_node in keys(placement_restriction)
            @constraint(model, sum(x[v_network, v_node, s_node] for s_node in placement_restriction[v_node]) == 1)
        end
    end


    ### =========== Nodes constraints

    # one substrate node per virtual node
    for v_network in instance.v_networks
        for v_node in vertices(v_network)
            @constraint(model, sum(x[v_network, v_node, s_node] for s_node in vertices(instance.s_network)) == 1)
        end
    end

    # if one to one : one virtual node per substrate node
    if one_to_one
        for s_node in vertices(instance.s_network)
            for v_network in instance.v_networks
                @constraint(model, sum(x[v_network, v_node, s_node] for v_node in vertices(v_network)) <= 1)
            end
        end
    end

    # node capacity
    for s_node in vertices(instance.s_network)
        @constraint(model, 
            sum( v_network[v_node][:dem] * x[v_network, v_node, s_node] 
                for v_network in instance.v_networks for v_node in vertices(v_network) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end


    ### ========== Edges constraints 
    
    # edge capacity
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge] 
                for v_network in instance.v_networks for v_edge in edges(v_network)) 
            <= 
            instance.s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(instance.s_network)
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                @constraint(model, 
                    x[v_network, src(v_edge), s_node] - x[v_network, dst(v_edge), s_node] 
                    <=
                    sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) - 
                        sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node))
                )
            end
        end
    end


    ### ========= Additional constraints
    
    ## Departure constraint
    if one_to_one
        if departure_cst
            for s_node in vertices(instance.s_network)
                for v_network in instance.v_networks
                    for v_node in vertices(v_network)
                        for v_edge in get_out_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                        end
                        for v_edge in get_in_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                        end
                    end
                end
            end
        end
    end

    ## Symmetric edges routings (because undirected !)
    undirected = true
    if undirected
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                for s_edge in edges(instance.s_network)
                    @constraint(model, y[v_network, v_edge, s_edge] == y[v_network, get_edge(v_network, dst(v_edge), src(v_edge)), 
                                                get_edge(instance.s_network, dst(s_edge), src(s_edge))])
                end
            end
        end
    end

    # Solving
    set_silent(model)
    optimize!(model)

    status = termination_status(model)

    if status != MOI.OPTIMAL
        println("Infeasible or unfinished: $status")
        return
    end

    print("Bah alors ?")
    # Get the solution
    x_values = value.(model[:x])
    y_values = value.(model[:y])
    columns = []
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            for s_node in vertices(instance.s_network)
                if x_values[v_network, v_node, s_node] > 0.99
                    append!(node_placement, s_node)
                end
            end
        end

        edge_routing = Dict()
        for v_edge in edges(v_network)
            if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
                edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
            end
            used_edges = []
            for s_edge in edges(instance.s_network)
                if y_values[v_network, v_edge, s_edge] > 0.99
                    push!(used_edges, s_edge)
                end
            end
            edge_routing[v_edge] = order_path(instance.s_network, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
        end
        m = Mapping(v_network, instance.s_network, node_placement, edge_routing)

        println(m)

        column = Column(m, objective_value(model))
        push!(columns, column)
    end

    return columns[1]
end






### ========= Column Generation result

# !!!!!!!!! This is not ready ! it's not taking into account a lot of things...



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
        end
        println("   For subgraph $(subgraph.graph), cost: $(current_subgraph_costs)")
        overall_subgraph_costs += current_subgraph_costs
    end

    println("Overall mapping of subgraphs costs : $(overall_subgraph_costs)")


    println("Master node placement:")
    for (v_node, placement) in mapping.master_placement
        println("   For vnode $(v_node)")
        for (s_node, val) in placement
            println("           $(s_node)  : $(val)")
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
end

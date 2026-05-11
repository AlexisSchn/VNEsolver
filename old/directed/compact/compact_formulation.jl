using Revise, JuMP, CPLEX, Gurobi
includet("../../../utils/import_utils.jl")
includet("solution_compact_fractional.jl")



struct Compact_Formulation
    model
    x_variables
    y_variables
end


function set_up_compact_model_gurobi(instance, one_to_one = false, departure_cst = false, symmetric = false)
    print("Constructing compact model... ")



    #### Model
    model = Model(Gurobi.Optimizer)
    #set_optimizer_attribute(model, "VarBranch", 3)
    #set_optimizer_attribute(model, "DisplayInterval", 1)
    #set_optimizer_attribute(model, "Seed", 0)

    ### Variables
    x_variables = @variable(model, x[v_network in instance.v_networks, vertices(v_network), vertices(instance.s_network)], binary=true);
    y_variables = @variable(model, y[v_network in instance.v_networks, edges(v_network), edges(instance.s_network)], binary=true);
    
    
    # BRANCHING PRIORITY

    for v_network in instance.v_networks
        for v_edge in edges(v_network)
            for s_edge in edges(instance.s_network)
                MOI.set(model, Gurobi.VariableAttribute("BranchPriority"), y[v_network, v_edge, s_edge], 0)
            end
        end
    end 

    for v_network in instance.v_networks
        for v_node in vertices(v_network)
            for s_node in vertices(instance.s_network)
                MOI.set(model, Gurobi.VariableAttribute("BranchPriority"), x[v_network, v_node, s_node], 1)
            end
        end
    end 



    
    #= TEST EARLY BRANCHING
    partition_sn = [    [8, 5, 9, 6, 15, 18, 17, 16, 10, 7, 2, 20],
                        [6, 14, 18, 19, 12, 21, 13, 1, 4],
                        [2, 20, 11, 19, 18, 15, 6, 9, 5, 8],
                        [10, 8, 7, 2, 20, 13, 21, 12, 19],
                        [18, 15, 17, 16, 14, 13, 1, 4, 3],
                        [8, 14, 21, 19, 9, 3]
    ]
    
    
    partition_sn = [    [8, 5, 9, 6, 15, 18, 17, 16, 10, 7, 2, 20],
                        [6, 14, 18, 19, 12, 21, 13, 1, 4]]

    partition_sn = [    [11, 19, 12, 21, 13, 14, 3, 4, 1],
                        [7, 2, 20, 11, 12, 21, 14, 13, 1, 3, 4],
                        [5, 8, 10, 7, 2, 20, 14, 3, 4, 1, 13]]


    v_node_central = [2, 11, 9, 7, 17, 13]

    p_variables = @variable(model, p[v_network in instance.v_networks, v_node in v_node_central, part in partition_sn], binary=true)

    for v_network in instance.v_networks
        for v_node in v_node_central
            for part in partition_sn
                MOI.set(model, Gurobi.VariableAttribute("BranchPriority"), p[v_network, v_node, part], 100)
            end
        end
    end 


    ############ LA SHEESH CONTRAINTE SUPPLEMENTAIRE
    for v_network in instance.v_networks
        for v_node in v_node_central
            for part in partition_sn
                print(part)
                @constraint(model, sum(x[v_network, v_node, s_node] for s_node in part) == p[v_network, v_node, part])
            end
        end
    end

    =#
    

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_network, v_node, s_node] 
        for v_network in instance.v_networks for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( instance.s_network[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge]
        for v_network in instance.v_networks for v_edge in edges(v_network) for s_edge in edges(instance.s_network) ))
    @objective(model, Min, placement_cost + routing_cost);


    ### Constraints
    


    ## Nodes

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


    ## Edges 
    
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
                    ==
                    sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) - 
                        sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node))
                )
            end
        end
    end


    ## Additional constraints : Node + Edge
    if one_to_one
        if departure_cst
            for s_node in vertices(instance.s_network)
                for v_network in instance.v_networks
                    for v_node in vertices(v_network)
                        for v_edge in get_out_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node)) <= 1 - x[v_network, v_node, s_node])
                        end
                        for v_edge in get_in_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                        end


                    end
                end
            end
        end
    end

    ## Symmetric edges (mostly for undirected)
    if symmetric
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                for s_edge in edges(instance.s_network)
                    @constraint(model, y[v_network, v_edge, s_edge] == y[v_network, get_edge(v_network, dst(v_edge), src(v_edge)), 
                                                get_edge(instance.s_network, dst(s_edge), src(s_edge))])
                end
            end
        end
    end
    
    println("done.")

    return Compact_Formulation(model, x_variables, y_variables)
end



function set_up_compact_model(instance, one_to_one = false, departure_cst = false, symmetric = false, silent = true)
    
    if !silent
        print("Constructing compact model... ")
    end
    #### Model
    model = Model(CPLEX.Optimizer)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Strategy_VariableSelect", 0)
    set_optimizer_attribute(model, "CPXPARAM_MIP_Display", 2)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Interval", 50)

    # CUTS
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Cuts_Nodecuts", 0)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Cuts_ZeroHalfCut", -1)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Cuts_LiftProj", -1)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Cuts_Gomory", -1)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Cuts_GUBCovers", 2)

    ### Variables
    x_variables = @variable(model, x[v_network in instance.v_networks, vertices(v_network), vertices(instance.s_network)], binary=true);
    y_variables = @variable(model, y[v_network in instance.v_networks, edges(v_network), edges(instance.s_network)], binary=true);
    #x_variables = @variable(model, 0 <=  x[v_network in instance.v_networks, vertices(v_network), vertices(instance.s_network)] <= 1);
    #y_variables = @variable(model, 0 <= y[v_network in instance.v_networks, edges(v_network), edges(instance.s_network)] <= 1);

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_network, v_node, s_node] 
        for v_network in instance.v_networks for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( instance.s_network[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge]
        for v_network in instance.v_networks for v_edge in edges(v_network) for s_edge in edges(instance.s_network) ))
    @objective(model, Min, placement_cost + routing_cost);


    ### Constraints

    ## Nodes

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


    ## Edges 
    
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


    ## Additional constraints : Node + Edge
    if one_to_one
        if departure_cst
            for s_node in vertices(instance.s_network)
                for v_network in instance.v_networks
                    for v_node in vertices(v_network)
                        for v_edge in get_out_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                            #@constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) <= x[v_network, v_node, s_node] + 0.01)
                        end
                        for v_edge in get_in_edges(v_network, v_node)
                            #@constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                        end
                    end
                end
            end
        end
    end

    ## Symmetric edges (mostly for undirected)
    if symmetric
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                for s_edge in edges(instance.s_network)
                    @constraint(model, y[v_network, v_edge, s_edge] == y[v_network, get_edge(v_network, dst(v_edge), src(v_edge)), 
                                                get_edge(instance.s_network, dst(s_edge), src(s_edge))])
                end
            end
        end
    end
    
    if !silent
        println("done.")
    end
    return Compact_Formulation(model, x_variables, y_variables)
end


function get_solution(instance, x_values, y_values)

    mappings = []
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
        push!(mappings, m)
    end

    return mappings
end


function solve_directed_compact_integer(instance, one_to_one = false, departure_cst = false, time_solver = 30, silent = true)

    # Set up the problem
    problem = set_up_compact_model(instance, one_to_one, departure_cst)


    # Solving
    print("Starting solving... ")
    set_time_limit_sec(problem.model, time_solver)
    if silent
        set_silent(problem.model)
    end
    optimize!(problem.model)
    println("done. Solving state: " * string(termination_status(problem.model)) * ", obj value: " * string(objective_value(problem.model)) * ", bound value: " * string(objective_bound(problem.model)))

    # Get the solution
    x_values = value.(problem.model[:x])
    y_values = value.(problem.model[:y])
    mappings = get_solution(instance, x_values, y_values)
    
    return mappings
end






function solve_directed_compact_fractional(instance, one_to_one = false, departure_cst = false, time_solver = 30, silent = true)

    # Set up the problem
    problem = set_up_compact_model(instance, one_to_one, departure_cst)
    relax_integrality(problem.model)

    #@constraint(problem.model, sum(problem.model[:y][v_network, v_edge, s_edge] for v_network in instance.v_networks for v_edge in edges(v_network) for s_edge in edges(instance.s_network)) >= 80 )
    model = problem.model
    #@constraint(model, + model[:x][instance.v_networks[1], 1, 2] - model[:x][instance.v_networks[1], 2, 2] + model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 2, 3)]  == 0)
    #@constraint(model, - model[:x][instance.v_networks[1], 1, 1] + model[:x][instance.v_networks[1], 2, 1] + model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 3, 1)]  == 0)
    #@constraint(model, - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 2)]  <= -1)
    #@constraint(model, - model[:x][instance.v_networks[1], 1, 2] + model[:x][instance.v_networks[1], 2, 2] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 2, 3)]  == -1)
    #@constraint(model, + model[:x][instance.v_networks[1], 1, 1] - model[:x][instance.v_networks[1], 2, 1] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 3, 1)]  == -1)
    #@constraint(model, + model[:x][instance.v_networks[1], 2, 1] + model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)]  ≤ 1)
    #@constraint(model, + model[:x][instance.v_networks[1], 1, 2] + model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)]  ≤ 1)
    #@constraint(model, - model[:x][instance.v_networks[1], 1, 2] - model[:x][instance.v_networks[1], 2, 1] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)]  ≤ -1)
    #@constraint(model, - model[:x][instance.v_networks[1], 1, 1] - model[:x][instance.v_networks[1], 2, 2] + model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)]  ≤ 0)    # Solving
    if silent
        set_silent(problem.model)
    end
    set_time_limit_sec(problem.model, time_solver)
    print("Starting solving... ")
    optimize!(problem.model)

    # Get the solution
    x_values = value.(problem.x_variables)
    y_values = value.(problem.y_variables)
    mappings = []
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            push!(node_placement, [])
            for s_node in vertices(instance.s_network)
                push!(node_placement[v_node], x_values[v_network, v_node, s_node])
            end
        end
        edge_routing = Dict()
        total_val = 0
        for v_edge in edges(v_network)
            edge_routing[v_edge] = Dict()
            val = 0
            for s_edge in edges(instance.s_network)
                val += y_values[v_network, v_edge, s_edge]
                edge_routing[v_edge][s_edge] = y_values[v_network, v_edge, s_edge]
            end
            total_val += val
            println("For edge $v_edge we have $val")
        end
        println("Total val : $total_val")
        m = MappingCompactFractional(v_network, instance.s_network, node_placement, edge_routing)
        push!(mappings, m)
    end

    return mappings
end





# Ne marche qu'en one to one pour l'instant...
function solve_directed_compact_with_subgraphs_constraints(instance, node_partitionning, time_solver = 30)
    
    # Set up the problem
    problem = set_up_compact_model(instance, true, true, true)

    
    # Subgraph stuff....
    for (i_vn, vn) in enumerate(instance.v_networks)
        println("For vn$i_vn...")
        for (i_subgraph, v_nodes) in enumerate(node_partitionning[i_vn])
            print("For subgraph $i_subgraph... ")
            subgraph = my_induced_subgraph(vn, v_nodes, "subgraph_$i_subgraph")
            # Résolution
            subinstance = InstanceVNE([subgraph], instance.s_network)
            mapping = solve_directed_compact_integer(subinstance, true, true, 60, false)[1]
            # get number of edges
            min_routing = 0
            for v_edge in edges(subgraph)
                min_routing += length(mapping.edge_routing[v_edge].edges)
            end
            print("Min edges: $min_routing instead of $(length(edges(subgraph)))\n")
    
            # adding the constraints
            @constraint(problem.model,
                sum( problem.model[:y][vn, get_edge(vn, v_nodes[src(v_edge)], v_nodes[dst(v_edge)]), s_edge] for v_edge in edges(subgraph) for s_edge in edges(instance.s_network))
                >=
                min_routing)
            
        end

    end 
    




    # Solving
    print("Starting solving... ")
    set_time_limit_sec(problem.model, time_solver)
    optimize!(problem.model)
    println("done. Solving state: " * string(termination_status(problem.model)) * ", obj value: " * string(objective_value(problem.model)) * ", bound value: " * string(objective_bound(problem.model)))

    # Get the solution
    x_values = value.(problem.model[:x])
    y_values = value.(problem.model[:y])
    mappings = get_solution(instance, x_values, y_values)
    
    return mappings

end

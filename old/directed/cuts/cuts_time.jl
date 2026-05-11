
using Revise, JuMP, CPLEX, Gurobi
includet("../../../utils/import_utils.jl")
includet("../compact/compact_formulation.jl")




function Base.show(io::IO, subcycle::SubCycle)

    println(io,
    "Subcycle $(subcycle.graph[]) with:
            graph $(subcycle.graph)
            assignment_vn_to_cycle: $(subcycle.assignment_vn_to_cycle)
            edges of original graph: $(subcycle.edges_of_original_graph)")

end


function set_up_compact_model_cuts(instance, one_to_one = true, departure_cst = true, symmetric = true)
    print("Constructing compact model... ")

    #### Model
    model = Model(CPLEX.Optimizer)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Strategy_VariableSelect", 0)
    set_optimizer_attribute(model, "CPXPARAM_MIP_Display", 2)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Interval", 50)

    # 
    set_optimizer_attribute(model, "CPXPARAM_MIP_Cuts_Nodecuts", -1)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Cuts_ZeroHalfCut", -1)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Cuts_LiftProj", -1)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Cuts_Gomory", -1)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Cuts_GUBCovers", 2)

    #set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Presolve", false)

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




function create_callback_cycles(problem, instance)


    return function my_callback_function(cb_data)

        #println("Starting callback...")


        x_values = callback_value.(cb_data, problem.model[:x])
        y_values = callback_value.(cb_data, problem.model[:y])

        println("I got the values !")
        for (i_vn, v_network) in enumerate(instance.v_networks)

            time1 = time()
            partial_placement = Dict()
            node_placement = []
            for v_node in vertices(v_network)
                push!(node_placement, [])
                for s_node in vertices(instance.s_network)
                    val = x_values[v_network, v_node, s_node]
                    if val > 0.9999
                        #println("Wow look ! $v_node is placed on $s_node")
                        partial_placement[v_node] = s_node
                    end
                    push!(node_placement[v_node], val)
                end
            end

            
            edge_routing = Dict()
            for v_edge in edges(v_network)
                edge_routing[v_edge] = Dict()
                for s_edge in edges(instance.s_network)
                    edge_routing[v_edge][s_edge] = y_values[v_network, v_edge, s_edge]
                end
            end
            time2 = time()
            time_spent_calling += (time2-time1)
            
            #m = MappingCompactFractional(v_network, instance.s_network, node_placement, edge_routing)

            #println("Current frac value: $(m.node_placement_cost + m.edge_routing_cost)")

            time3 = time()
            time_spent_solving += (time3 - time2)
            println("Time spent calling : $time_spent_calling, time spent solving: $time_spent_solving")
            
            # verify that it does improve 
        end   
         
    end
end

function solve_directed_cuts(instance, one_to_one = false, departure_cst = false, time_solver = 60, silent = false)

    # Set up the problem
    problem = set_up_compact_model_cuts(instance, one_to_one, departure_cst)

    #println("CYCLES: $cycles")
    # Solving
    print("Starting solving... ")
    set_time_limit_sec(problem.model, 500)
    if silent
        set_silent(problem.model)
    end
    my_callback_function = create_callback_cycles(problem, instance)
    set_optimizer_attribute(problem.model, "CPXPARAM_Threads", 1)
    set_attribute(problem.model, MOI.UserCutCallback(), my_callback_function)
    
    optimize!(problem.model)
    println("done. Solving state: " * string(termination_status(problem.model)) * ", obj value: " * string(objective_value(problem.model)) * ", bound value: " * string(objective_bound(problem.model)))

    # Get the solution
    #x_values = value.(problem.model[:x])
    #y_values = value.(problem.model[:y])
    #mappings = get_solution(instance, x_values, y_values)
    
    return mappings
end



# 2 : make a solver, with fixed nodes

# 3 : observe, perhaps add the constraint ?

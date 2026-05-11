#using Revise, JuMP, CPLEX, Gurobi
using Revise, JuMP, CPLEX
using GraphsFlows

includet("../utils/import_utils.jl")




# ========== CLASSICAL STUFF

function set_up_problem_ff_extended(instance, model)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    ### Variables
    @variable(model, x[vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);
    @variable(model, 1 >= z[edges(v_network), vertices(s_network), vertices(s_network)] >= 0)

    

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);


    ###=========== Constraints

    ##---- Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(instance.s_network)) == 1)
    end

    # one to one : one virtual node per substrate node
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end

    # node capacity : NOT USELESS AHHHHHHHHh
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(v_network[v_node][:dem] * x[v_node, s_node] for v_node in vertices(v_network)) <= sum(s_network[s_node][:cap]))
    end


    ##--- Edges 
    
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

    
    ## Departure constraints    
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end
    
    # Outgoing edges cap: pretty stupid but useful
    for v_node in vertices(v_network)
        for s_node in vertices(s_network)
            v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
            necessary_bw = 0 + sum(v_network[src(v_edge), dst(v_edge)][:dem] for v_edge in v_edges_incident)

            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = 0 +sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end
    

    # zizi Variables

    # zizi must correspond to the x variables
    for v_edge in edges(v_network)
        for s_node_1 in vertices(s_network)
            @constraint(model, sum( z[v_edge, s_node_1, s_node_2] for s_node_2 in vertices(s_network) ) >= x[src(v_edge), s_node_1])
            @constraint(model, z[v_edge, s_node_1, s_node_1] == 0)
        end
    end
    for v_edge in edges(v_network)
        for s_node_2 in vertices(s_network)
            @constraint(model, sum( z[v_edge, s_node_1, s_node_2] for s_node_1 in vertices(s_network) ) >= x[dst(v_edge), s_node_2])
            @constraint(model, z[v_edge, s_node_2, s_node_2] == 0)
        end
    end

    # one zizi per v edge: nope


end



function solve_compact_extended(instance; time_solver = 30, stay_silent=false, linear=false)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir


    model = Model(CPLEX.Optimizer)
    set_up_problem_ff_extended(instance, model)

    set_time_limit_sec(model, time_solver)
    if stay_silent
        set_silent(model)
    else
        print("Starting solving... ")
    end


    my_callback_function = create_callback_connective_cuts(model, instance)
    set_optimizer_attribute(model, "CPXPARAM_Threads", 1)

    set_attribute(model, MOI.UserCutCallback(), my_callback_function)
    
    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Infeasible or unfinished: $status")
        return -999, objective_bound(model), node_count(model)
    end

    println("nb de noeuds: $(node_count(model))")
    println("Lower bound: $(objective_bound(model))")
    println("The objective is : $(objective_value(model))")

    #=
    if !stay_silent

        x_values = value.(model[:x])
        y_values = value.(model[:y])
    
        println("Node placement:")
        for v_node in vertices(v_network)
            for s_node in vertices(s_network_dir)
                if x_values[v_node, s_node] > 0.01
                    println("$v_node is placed on $s_node")
                end
            end
        end
        println("\nEdge routing:")
        for v_edge in edges(v_network)
            print("Routing of $v_edge : ")
            for s_edge in edges(s_network_dir)
                if y_values[v_edge, s_edge] > 0.01
                    print(" $s_edge")
                end
            end
            print("\n")
        end
    end
    =#

    result = objective_value(model)
    lb = objective_bound(model)
    nbnodes = node_count(model)
    return result, lb, nbnodes
end




function solve_compact_extended_chilling(instance; time_solver = 30, stay_silent=false, linear=false)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir


    model = Model(CPLEX.Optimizer)
    set_up_problem_ff_extended(instance, model)

    set_time_limit_sec(model, time_solver)
    if stay_silent
        set_silent(model)
    else
        print("Starting solving... ")
    end


    #my_callback_function = create_callback_connective_cuts(model, instance)
    #set_attribute(model, MOI.UserCutCallback(), my_callback_function)

    set_optimizer_attribute(model, "CPXPARAM_Threads", 1)

    
    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Infeasible or unfinished: $status")
        return -999, objective_bound(model), node_count(model)
    end

    println("nb de noeuds: $(node_count(model))")
    println("Lower bound: $(objective_bound(model))")
    println("The objective is : $(objective_value(model))")

    #=
    if !stay_silent

        x_values = value.(model[:x])
        y_values = value.(model[:y])
    
        println("Node placement:")
        for v_node in vertices(v_network)
            for s_node in vertices(s_network_dir)
                if x_values[v_node, s_node] > 0.01
                    println("$v_node is placed on $s_node")
                end
            end
        end
        println("\nEdge routing:")
        for v_edge in edges(v_network)
            print("Routing of $v_edge : ")
            for s_edge in edges(s_network_dir)
                if y_values[v_edge, s_edge] > 0.01
                    print(" $s_edge")
                end
            end
            print("\n")
        end
    end
    =#

    result = objective_value(model)
    lb = objective_bound(model)
    nbnodes = node_count(model)
    return result, lb, nbnodes
end





function create_callback_connective_cuts(model, instance)


    i = 0

    return function my_callback_function(cb_data)

        if i==0 
            v_network = instance.v_network
            s_network = instance.s_network
            s_network_dir = instance.s_network_dir

            println("Starting callback...")

            x_values = callback_value.(cb_data, model[:x])
            y_values = callback_value.(cb_data, model[:y])
            nb_cons = 0

            #println("Edges: $(edges(v_network))")
            for v_edge in edges(v_network)
                #println("Doing edge $v_edge")
                
                weights_from_pl = zeros(Float32, nv(s_network), nv(s_network))
                for s_edge in edges(instance.s_network_dir)
                    weights_from_pl[src(s_edge), dst(s_edge)] = abs(y_values[v_edge, s_edge])
                end

                #println("the sum is: $(sum(weights_from_pl))")

                #=
                
                println("Look at my beautiful matrix uwu:")
                for i in 1:nv(instance.s_network)
                    for j in 1:nv(instance.s_network)
                        print("$(weights_from_pl[(i-1)*(nv(instance.s_network))+j]) ")
                    end
                    println()
                end
                =#

                #println(weights_from_pl)
                
                #print(weights_from_pl)

                for s_node_1 in vertices(s_network)
                    x_val_1 = x_values[src(v_edge), s_node_1]
                    for s_node_2 in vertices(s_network)
                        
                        if s_node_1 != s_node_2
                            x_val_2 = x_values[dst(v_edge), s_node_2]

            
                            if x_val_1 > 0.03
                                if x_val_2 > 0.03
                                    # get the max flow on the correct graph...
                                    #trucmuche = mincut(s_network, s_node_1, s_node_2)
                                    #print("For $v_edge, from $s_node_1 to $s_node_2, , we should have $(min(x_val_1, x_val_2))") 
                                    part1, part2, value = GraphsFlows.mincut(instance.s_network_dir.graph, s_node_1, s_node_2, weights_from_pl, GraphsFlows.DinicAlgorithm())
                                    
                                    # get a min cut from it
                                    
                                    if value < min(x_val_1, x_val_2)-0.05
                                        #print(" but we got $value")
                                        #print(" And the parts are $part1, and $part2")
                                        cut_s_edges = get_edges_from_S1_to_S2(s_network_dir, part1, part2)
                                        #print(" and my cut edges: $cut_s_edges")
                                        @constraint(model, sum(model[:y][v_edge, s_edge] for s_edge in cut_s_edges)
                                            ≥ 
                                            model[:z][v_edge, s_node_1, s_node_2])

                                        nb_cons += 1
                                    end

                                end
                            end


                            # then add the constraint...
                        end

                    end
                end
            end

            println("Finished, $nb_cons cuts added")
        end

        i+=1
    end
end




function get_edges_from_S1_to_S2(graph, S1, S2)

    cut_edges = []
    for edge in edges(graph)
        if src(edge) ∈ S1 && dst(edge) ∈ S2
            push!(cut_edges, edge)
        end
    end

    return cut_edges
end





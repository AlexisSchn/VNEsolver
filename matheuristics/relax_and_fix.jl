
using Revise, JuMP, CPLEX, Gurobi
includet("../../../utils/import_utils.jl")




function set_up_problem_placement_restrict(instance, model, placement_restriction, nodes_to_restrict)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    ### Variables
    @variable(model, 0<=x[vertices(v_network), vertices(s_network)] <=1);
    for s_node in vertices(s_network)
        for v_node in nodes_to_restrict
            set_integer(x[v_node, s_node])
        end
    end
    @variable(model, 0 <= y[edges(v_network), edges(s_network_dir)] <= 1);

    

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
        if v_node ∈ keys(placement_restriction)   
            @constraint(model, sum(x[v_node, s_node] for s_node in placement_restriction[v_node]) == 1)
        end
    end

    # one to one : one virtual node per substrate node
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end

    # node capacity : NOT USELESS AHHHHHHHHh
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(v_network[v_node][:dem] * x[v_node, s_node] for v_node in vertices(v_network)) <= sum(s_network[s_node][:cap]))
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

    #--------- departure constraint
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end
    
    #--------------- infeasibility constraints: if there is just not enough capacity on edges yknow...
    nb_var_less = 0
    for v_node in vertices(v_network)
        for s_node in vertices(s_network)
            v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
            necessary_bw = sum(v_network[src(v_edge), dst(v_edge)][:dem] for v_edge in v_edges_incident)

            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident)
            if necessary_bw > available_bw
                nb_var_less += 1
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end

        
end



function set_up_problem_placement_restrict_end(instance, model, placement_restriction)

    v_network = instance.v_network
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
        if v_node ∈ keys(placement_restriction)   
            @constraint(model, sum(x[v_node, s_node] for s_node in placement_restriction[v_node]) == 1)
        end
    end

    # one to one : one virtual node per substrate node
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end

    # node capacity : NOT USELESS AHHHHHHHHh
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(v_network[v_node][:dem] * x[v_node, s_node] for v_node in vertices(v_network)) <= sum(s_network[s_node][:cap]))
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

    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end
    
end




function relax_and_fix(instance)

    v_network = instance.v_network
    s_network = instance.s_network

    time_start = time()

    nb_neighbors = []
    nb_neighbors_done = []
    for v_node in vertices(v_network)
        push!(nb_neighbors, length(neighbors(v_network, v_node)))
        push!(nb_neighbors_done, 0)
    end


    restrictions = Dict()
    restricted_nodes = []

    # we simply takes the most dense node and make it start here
    sorted_degree = sort(collect(vertices(v_network)), by=v_node -> -degree(v_network, v_node))
    nodes_to_restrict_next = sorted_degree[1:1]

    iter=0
    while length(restricted_nodes) != length(vertices(v_network))
        iter+=1
        current_model = Model(CPLEX.Optimizer)
        println("\n---------------\nIter $iter\n---------------")
        #println("Current placement restriction: $(restrictions)")
        set_up_problem_placement_restrict(instance, current_model, restrictions, nodes_to_restrict_next) 
        #add_star_constraints(instance, current_model)
        set_silent(current_model)
    
        optimize!(current_model)

        if !has_values(current_model)
            println("Infeasible model !")
            return
        end
        #println("Solving time: $(solve_time(current_model)) ")
        println("Value LP: $(objective_value(current_model)), solved in $(solve_time(current_model))")
        

        x_values = value.(current_model[:x])

        print("New placement:")
        add_all_vnodes = true
        for s_node in vertices(s_network)
            for v_node in vertices(v_network)
                if x_values[v_node, s_node] > 0.99
                    if v_node in nodes_to_restrict_next
                        restrictions[v_node] = s_node
                        push!(restricted_nodes, v_node)
                        print(" $v_node => $s_node")
                    elseif add_all_vnodes && (v_node ∉ restricted_nodes)
                        restrictions[v_node] = s_node
                        push!(restricted_nodes, v_node)
                        print(" $v_node => $s_node")
                    end
                end
            end
        end
        print("\n")

        # -------- NEXT NODES SELECTION: quite tricky ?
        # idea : 1) most central left, 2) most neighbors 3) less frac 4) all neighbors ?
        next_nodes_choice = 2
        if length(restricted_nodes) != length(vertices(v_network))

            if next_nodes_choice == 1
                ratio = Dict()
                for v_node in vertices(v_network)
                    if v_node ∉ restricted_nodes
                        #ratio[v_node] = -nb_neighbors_done[v_node] / nb_neighbors[v_node] # with less neighbor already fixed first
                        #ratio[v_node] = nb_neighbors_done[v_node] / nb_neighbors[v_node] # with most neighbor already fixed first
                        ratio[v_node] = degree(v_network, v_node) # with highest degree node
                    end
                end
                sorted_keys = sort(collect(keys(ratio)), by=x -> -ratio[x])
                nodes_to_restrict_next = [sorted_keys[1]]
            elseif next_nodes_choice == 2

                nodes_to_restrict_next= []
                for v_node in restricted_nodes
                    for neighbor in neighbors(v_network, v_node)
                        if neighbor ∉ restricted_nodes
                            if neighbor ∉ nodes_to_restrict_next
                                push!(nodes_to_restrict_next, neighbor)
                            end
                        end
                    end
                end
            end

            println("Next nodes to restrict: $nodes_to_restrict_next")


        end
           

    end

    #total_time = time() - time0
    #println("Total time: $(total_time), in lp : $(total_solve_time), const: $(time_constructing_model), getting sols etc: $(time_choosing_x)")

    
    # Final resolution
    println("\n---------------\nFinal iteration\n---------------")
    println("Placement: $restrictions")
    model_end = Model(CPLEX.Optimizer)
    set_up_problem_placement_restrict_end(instance, model_end, restrictions)
    set_silent(model_end)
    optimize!(model_end)
    println("ILP: $(objective_value(model_end)), solved in $(solve_time(model_end))\n")



    time_algo = time() - time_start
    println("Overall, it took $time_algo to solve the instance, for a solution of $(objective_value(model_end))")

end


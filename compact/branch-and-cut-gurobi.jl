
#using Revise, JuMP, CPLEX, Gurobi
using Revise, JuMP, Gurobi
using Graphs, GraphsFlows

includet("../utils/import_utils.jl")




# ========== CLASSICAL STUFF
function set_up_problem_ff_plus(instance, model)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    ### Variables
    @variable(model, x[vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);

    

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * y[v_edge, s_edge]
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
        @constraint(model, sum( x[v_node, s_node] for v_node in vertices(v_network)) <= sum(s_network[s_node][:cap]))
    end


    ##--- Edges 
    
    # edge capacity (undirected version !)
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
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
    
    
    
    #= Star capacity constraint => NOT FOR NOW?
    for v_node in vertices(v_network)
        necessary_bw = degree(v_network, v_node) 
        
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident;init=0.)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end
    =#  
end



function solve_branch_and_cut(instance; time_solver = 100)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir

    time_start = time()
    model = Model(Gurobi.Optimizer)
    set_up_problem_ff_plus(instance, model)

    set_time_limit_sec(model, time_solver)

    #set_optimizmer_attribute(model, "CPXPARAM_MIP_Strategy_Search", 1)
    #set_optimizer_attribute(model, "CPX_PARAM_REPEATPRESOLVE", 0)
    set_optimizer_attribute(model, "Threads", 1)

    # setting things for the branch and cut
    n_s = nv(s_network_dir)
    augmented_network = copy(s_network_dir.graph)
    add_vertex!(augmented_network) # n_s + 1: t(\ebar)
    for i in 1:n_s
        add_edge!(augmented_network, i, n_s+1)
    end

    nb_cuts_overall = 0
    time_seperation_overall = 0
    function find_cutset_cuts(cb_data)
        
        time_start_separation = time()
        nb_cuts_current = 0

        x_values = callback_value.(cb_data, model[:x])
        y_values = callback_value.(cb_data, model[:y])


        for v_edge in edges(v_network)
            matrix_flows = zeros(n_s+1, n_s+1)
 
            for s_edge in edges(s_network_dir)
                matrix_flows[src(s_edge), dst(s_edge)] = y_values[v_edge, s_edge]
            end

            for s_node in vertices(s_network_dir)
                x_values[src(v_edge), s_node] < 0.001 && continue
                # Finish flow matrix for that node - remember that no flow between u and t(\ebar) here
                for other_s_node in vertices(s_network_dir)
                    if s_node == other_s_node
                        matrix_flows[other_s_node, n_s+1] = 0.
                    else 
                        matrix_flows[other_s_node, n_s+1] = x_values[dst(v_edge), other_s_node]
                    end
                end
                
                #println("Matrix flow : $matrix_flows")
                # min cut between s_node and n_s+1
                (part1, part2, flow) = GraphsFlows.mincut(augmented_network, s_node, n_s+1, matrix_flows, DinicAlgorithm())

                if flow < x_values[src(v_edge), s_node] - 0.0001 
                    cut_s_edges = get_edges_from_S1_to_S2(s_network_dir, part1, part2)
                    cut = @build_constraint(sum(model[:y][v_edge, s_edge] for s_edge in cut_s_edges) + sum( model[:x][dst(v_edge), s_other_node] for s_other_node in part1)
                        >=
                        model[:x][src(v_edge), s_node] + model[:x][dst(v_edge), s_node]
                    )
                    MOI.submit(model, MOI.UserCut(cb_data), cut)

                    nb_cuts_current += 1
                    nb_cuts_overall += 1
                end
            end


        end

        time_seperation_overall += time() - time_start_separation
        #println("Added $nb_cuts_current this time.")
        #println("Added $nb_cuts_overall so far, in $time_overall.")

        return

    end


    MOI.set(model, MOI.UserCutCallback(), find_cutset_cuts) 

    
    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Infeasible or unfinished: $status")
        return ( sol_value= -1,
                lower_bound = -1.,
                gap = -1.,
                node_count = node_count(model),
                time_solving = (time() - time_start),
                time_seperation = time_seperation_overall,
                nb_cuts = nb_cuts_overall
        )
    end




    # Get the solution

    
    

    println("Find the solution $(objective_value(model)) in $(time() - time_start)")
    println("Stats branch and cut: $time_seperation_overall s to find $nb_cuts_overall cuts")


    return    ( sol_value= objective_value(model),
                lower_bound = objective_bound(model),
                gap = relative_gap(model),
                node_count = node_count(model),
                time_solving = (time() - time_start),
                time_seperation = time_seperation_overall,
                nb_cuts = nb_cuts_overall
    )
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








function solve_simple_gurobi(instance; time_solver = 100)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir

    time_start = time()
    model = Model(Gurobi.Optimizer)
    set_up_problem_ff_plus(instance, model)

    set_time_limit_sec(model, time_solver)

    # On
    set_optimizer_attribute(model, "Threads", 1)

    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Infeasible or unfinished: $status")
        return ( sol_value= -1,
            lower_bound = -1.,
            gap = -1.,
            node_count = node_count(model),
            time_solving = (time() - time_start)
        )
    end




    # Get the solution


    return ( sol_value= objective_value(model),
        lower_bound = objective_bound(model),
        gap = relative_gap(model),
        node_count = node_count(model),
        time_solving = (time() - time_start)
    )
end


#using Revise, JuMP, CPLEX, Gurobi
using Revise, JuMP, CPLEX
using Graphs, GraphsFlows
using SparseArrays


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
    
    
end


function set_up_problem_ff_plusplus(instance, model, continuity_degree, continuity_cap)

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

    nb_continuity = 0
    for s_edge_in in edges(s_network_dir)
        for v_edge in edges(v_network)
            s_node = dst(s_edge_in)
            
            if degree(s_network, s_node) < continuity_degree
                if continuity_cap  && (s_network[s_node][:cap]==0)
                    @constraint(model, sum(y[v_edge, s_edge_out] for s_edge_out in get_out_edges(s_network_dir, s_node) ) + x[dst(v_edge), s_node] 
                        >= y[v_edge, s_edge_in] + y[v_edge, get_reverse_edge(s_network_dir, s_edge_in)] )
                    nb_continuity += 1
                elseif !continuity_cap
                    @constraint(model, sum(y[v_edge, s_edge_out] for s_edge_out in get_out_edges(s_network_dir, s_node) ) + x[dst(v_edge), s_node] 
                        >= y[v_edge, s_edge_in] + y[v_edge, get_reverse_edge(s_network_dir, s_edge_in)] )
                    nb_continuity += 1
                end
            end
        end
    end
    
    #=
    # Star capacity constraint
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
        
    return nb_continuity
end




function solve_branch_and_cut(instance; solver="cplex", time_solver = 100, alpha_acceptance=0.05)
    
    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    time_start = time()
    if solver == "cplex"
        model = Model(CPLEX.Optimizer)
        set_up_problem_ff_plus(instance, model)
        set_optimizer_attribute(model, "CPXPARAM_RandomSeed", rand(1:10^9))  # or "CPX_PARAM_RANDOMSEED"
    elseif solver == "gurobi"
        model = Model(Gurobi.Optimizer)
        set_up_problem_ff_plus(instance, model)
        set_optimizer_attribute(model, "Threads", 1)
    else
        println("Solver didnt found. What?")
        return
    end

    set_time_limit_sec(model, time_solver)


    #set_optimizmer_attribute(model, "CPXPARAM_MIP_Strategy_Search", 1)
    #set_optimizer_attribute(model, "CPX_PARAM_REPEATPRESOLVE", 0)
    #
    #set_optimizer_attribute(model, "CPXPARAM_Parallel", -1)

    # setting things for the branch and cut
    n_s = nv(s_network_dir)
    augmented_network = copy(s_network_dir.graph)
    add_vertex!(augmented_network) # n_s + 1: t(\ebar)
    for i in 1:n_s
        if s_network[i][:cap] > 0
            add_edge!(augmented_network, i, n_s+1)
        end
    end



    x = model[:x]
    y = model[:y]

    nb_cuts_overall = 0
    nb_min_cut = 0
    nb_call = 0
    time_seperation_overall = 0
    time_min_cut = 0
    time_get_values = 0
    time_matrix_flow = 0
    time_constraint = 0
    
    matrix_flows = zeros(Float32, n_s+1, n_s+1)

    # Variable stuff
    n_x = nv(v_network) * nv(s_network)
    n_y = ne(v_network) *  ne(s_network_dir)
    x_idx = [MOI.VariableIndex(i) for i in 1:n_x]
    y_idx = [MOI.VariableIndex(i) for i in n_x+1:n_x+n_y]
    x_vals = zeros(Float64, n_x)    
    y_vals = zeros(Float64, n_y)
    backend = JuMP.backend(model)


    function find_cutset_cuts(cb_data)
        
        nb_call += 1
        if nb_call % 100 == 0
            println("At call $nb_call of callback, did $nb_min_cut mincuts and found $nb_cuts_overall cuts!")
        end


        if alpha_acceptance > 0.99 # Ahah don't worry (this is ugly but it's just for the tests)
            return
        end

        time_start_separation = time()
        nb_cuts_current = 0
        time_beg_getting_values = time()
        
        #=
        x_values = callback_value.(cb_data, x)
        y_values = callback_value.(cb_data, y)
        =#

        MOI.get!(x_vals,
            backend.optimizer.model,
            MOI.CallbackVariablePrimal(cb_data),
            x_idx
        )

        MOI.get!(y_vals,
            backend.optimizer.model,
            MOI.CallbackVariablePrimal(cb_data),
            y_idx
        )

        time_get_values += time()-time_beg_getting_values

        i_v_edge = 0
        for v_edge in edges(v_network)

            time_matrix_flow_beg = time()
            matrix_flows .= 0.0  # resets in place

            i_v_edge += 1
            i_s_edge = 0
            for s_edge in edges(s_network_dir)
                i_s_edge += 1
                matrix_flows[src(s_edge), dst(s_edge)] = y_vals[ ne(v_network) * (i_s_edge-1) + i_v_edge]
            end

            #@inbounds matrix_flows .= clamp.(matrix_flows, 0.0, 1.0)
            @inbounds @simd for i in eachindex(matrix_flows)
                v = matrix_flows[i]
                if !(isfinite(v)) || v < 1e-12
                    matrix_flows[i] = 0.0
                elseif v > 1.0
                    matrix_flows[i] = 1.0
                end
            end

            time_matrix_flow += time() - time_matrix_flow_beg

            #println("matrix flow: $matrix_flows")
            
            for s_node in vertices(s_network_dir)
                val = x_vals[nv(v_network) * (s_node-1) + src(v_edge)]
                
                val < alpha_acceptance * 1.01 && continue # reduces the number of separations
                val > 0.9 && continue # prevent bug in the mincut!

                nb_min_cut += 1

                # Finish flow matrix for that node - remember that no flow between u and t(\ebar) here
                for other_s_node in vertices(s_network_dir)
                    if s_node == other_s_node
                        matrix_flows[other_s_node, n_s+1] = 0.
                    else 
                        matrix_flows[other_s_node, n_s+1] = x_vals[nv(v_network) * (other_s_node-1) + dst(v_edge)]
                    end
                end
                
                # min cut between s_node and n_s+1
                time_beg_cut = time()
                (part1, part2, flow) = GraphsFlows.mincut(augmented_network, s_node, n_s+1, matrix_flows, DinicAlgorithm())
                time_min_cut += time()-time_beg_cut
                
                if flow < val - alpha_acceptance

                    if n_s+1 ∈ part1
                        println("WOWOWOW NS+1 IS IN PART1?? $part1 and $part2 wthhh")
                        println("I mean, I'm doing a cut between $s_node and $(n_s+1)... Is that algorithm stupid?")
                        break
                    end
                    if s_node ∈ part2
                        println("WOWOWOW $s_node IS IN PART2?? $part1 and $part2 wthhh")
                        println("I mean, I'm doing a cut between $s_node and $(n_s+1)... Is that algorithm stupid?")
                        break
                    end
                    
                    #=
                    println("Damn new cut for $v_edge on $s_node: $part1")
                    for s_node_other ∈ part1
                        if s_node_other != s_node
                            println("       with x_$(s_node_other): $(x_vals[nv(v_network) * (s_node_other-1) + dst(v_edge)])")
                        end
                    end
                    =#

                    
                    cut_s_edges = get_edges_from_S1_to_S2(s_network_dir, part1, part2) # Not the fastest, but not too bad.

                    #=
                    time_constraint_beg = time()

                    terms = Vector{MOI.ScalarAffineTerm{Float64}}()
                    # y part
                    for s_edge in cut_s_edges
                        push!(terms, MOI.ScalarAffineTerm(1.0, y[v_edge, s_edge].index))
                    end
                    # x[dst(v_edge), s_other_node] part
                    for s_other_node in part1
                        if s_other_node != s_node
                            push!(terms, MOI.ScalarAffineTerm(1.0, MOI.VariableIndex(nv(v_network) * (s_other_node-1) + dst(v_edge))))
                        end
                    end
                    push!(terms, MOI.ScalarAffineTerm(-1.0, MOI.VariableIndex(nv(v_network) * (s_node-1) + src(v_edge))))
                
                    
                
                    MOI.submit(backend, MOI.UserCut(cb_data),
                        MOI.ScalarAffineFunction(terms, 0.0),
                        MOI.GreaterThan(0.0)
                    )
                    

                    time_constraint += time() - time_constraint_beg

                    =#
            
                    
                    time_constraint_beg = time()

                    cut_s_edges = get_edges_from_S1_to_S2(s_network_dir, part1, part2) # Not the fastest, but not too bad.

                    cut = @build_constraint(sum(y[v_edge, s_edge] for s_edge in cut_s_edges) + sum( x[dst(v_edge), s_other_node] for s_other_node in part1)
                        >=
                        x[src(v_edge), s_node] + x[dst(v_edge), s_node]
                    ) 

                    MOI.submit(model, MOI.UserCut(cb_data), cut)
                    time_constraint += time() - time_constraint_beg


                    #println("Size of cut: $(length(part1))")
                    

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
                nb_cuts = nb_cuts_overall,
                nb_call = nb_call,
                nb_min_cut = nb_min_cut
        )
    end




    # Get the solution

    
    

    println("Find the solution $(objective_value(model)) in $(time() - time_start)")
    println("Stats branch and cut: $time_seperation_overall s overall, $time_min_cut on mincut, $time_get_values getting values, $time_matrix_flow on flow matrix, $time_constraint on constraint
                $nb_call calls, $nb_min_cut mincuts computed, $nb_cuts_overall cuts found")


    return    ( sol_value= objective_value(model),
                lower_bound = objective_bound(model),
                gap = relative_gap(model),
                node_count = node_count(model),
                time_solving = (time() - time_start),
                time_seperation = time_seperation_overall,
                nb_cuts = nb_cuts_overall, 
                nb_call = nb_call,
                nb_min_cut = nb_min_cut
    )
end



function solve_branch_and_cut_plusplus(instance; solver="cplex", time_solver = 100, alpha_acceptance=0.05, continuity_degree=1000, continuity_cap=true)
    
    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    time_start = time()
    model = Model(CPLEX.Optimizer)
    nb_continuity = set_up_problem_ff_plusplus(instance, model, continuity_degree, continuity_cap)
    set_optimizer_attribute(model, "CPXPARAM_RandomSeed", rand(1:10^9))  # or "CPX_PARAM_RANDOMSEED"

    set_time_limit_sec(model, time_solver)


    #set_optimizmer_attribute(model, "CPXPARAM_MIP_Strategy_Search", 1)
    #set_optimizer_attribute(model, "CPX_PARAM_REPEATPRESOLVE", 0)
    #
    #set_optimizer_attribute(model, "CPXPARAM_Parallel", -1)

    # setting things for the branch and cut
    n_s = nv(s_network_dir)
    augmented_network = copy(s_network_dir.graph)
    add_vertex!(augmented_network) # n_s + 1: t(\ebar)
    for i in 1:n_s
        if s_network[i][:cap] > 0
            add_edge!(augmented_network, i, n_s+1)
        end
    end



    x = model[:x]
    y = model[:y]

    nb_cuts_overall = 0
    nb_min_cut = 0
    nb_call = 0
    time_seperation_overall = 0
    time_min_cut = 0
    time_get_values = 0
    time_matrix_flow = 0
    time_constraint = 0
    
    matrix_flows = zeros(Float32, n_s+1, n_s+1)

    # Variable stuff
    n_x = nv(v_network) * nv(s_network)
    n_y = ne(v_network) *  ne(s_network_dir)
    x_idx = [MOI.VariableIndex(i) for i in 1:n_x]
    y_idx = [MOI.VariableIndex(i) for i in n_x+1:n_x+n_y]
    x_vals = zeros(Float64, n_x)    
    y_vals = zeros(Float64, n_y)
    backend = JuMP.backend(model)


    function find_cutset_cuts(cb_data)
        
        nb_call += 1
        if nb_call % 100 == 0
            println("At call $nb_call of callback, did $nb_min_cut mincuts and found $nb_cuts_overall cuts!")
        end


        if alpha_acceptance > 0.99 # Ahah don't worry (this is ugly but it's just for the tests)
            return
        end

        time_start_separation = time()
        nb_cuts_current = 0
        time_beg_getting_values = time()
        
        #=
        x_values = callback_value.(cb_data, x)
        y_values = callback_value.(cb_data, y)
        =#

        MOI.get!(x_vals,
            backend.optimizer.model,
            MOI.CallbackVariablePrimal(cb_data),
            x_idx
        )

        MOI.get!(y_vals,
            backend.optimizer.model,
            MOI.CallbackVariablePrimal(cb_data),
            y_idx
        )

        time_get_values += time()-time_beg_getting_values

        i_v_edge = 0
        for v_edge in edges(v_network)

            time_matrix_flow_beg = time()
            matrix_flows .= 0.0  # resets in place

            i_v_edge += 1
            i_s_edge = 0
            for s_edge in edges(s_network_dir)
                i_s_edge += 1
                matrix_flows[src(s_edge), dst(s_edge)] = y_vals[ ne(v_network) * (i_s_edge-1) + i_v_edge]
            end

            #@inbounds matrix_flows .= clamp.(matrix_flows, 0.0, 1.0)
            @inbounds @simd for i in eachindex(matrix_flows)
                v = matrix_flows[i]
                if !(isfinite(v)) || v < 1e-12
                    matrix_flows[i] = 0.0
                elseif v > 1.0
                    matrix_flows[i] = 1.0
                end
            end

            time_matrix_flow += time() - time_matrix_flow_beg

            #println("matrix flow: $matrix_flows")
            
            for s_node in vertices(s_network_dir)
                val = x_vals[nv(v_network) * (s_node-1) + src(v_edge)]
                
                val < alpha_acceptance * 1.01 && continue # reduces the number of separations
                val > 1-alpha_acceptance && continue # reduces the nuber of separatios, for free
                val > 0.99 && continue # prevent bug in the mincut!

                nb_min_cut += 1

                # Finish flow matrix for that node - remember that no flow between u and t(\ebar) here
                for other_s_node in vertices(s_network_dir)
                    if s_node == other_s_node
                        matrix_flows[other_s_node, n_s+1] = 0.
                    else 
                        matrix_flows[other_s_node, n_s+1] = x_vals[nv(v_network) * (other_s_node-1) + dst(v_edge)]
                    end
                end
                
                # min cut between s_node and n_s+1
                time_beg_cut = time()
                (part1, part2, flow) = GraphsFlows.mincut(augmented_network, s_node, n_s+1, matrix_flows, DinicAlgorithm())
                time_min_cut += time()-time_beg_cut
                
                if flow < val - alpha_acceptance

                    if n_s+1 ∈ part1
                        println("WOWOWOW NS+1 IS IN PART1?? $part1 and $part2 wthhh")
                        println("I mean, I'm doing a cut between $s_node and $(n_s+1)... Is that algorithm stupid?")
                        break
                    end
                    if s_node ∈ part2
                        println("WOWOWOW $s_node IS IN PART2?? $part1 and $part2 wthhh")
                        println("I mean, I'm doing a cut between $s_node and $(n_s+1)... Is that algorithm stupid?")
                        break
                    end
                    
                    
                    
                    cut_s_edges = get_edges_from_S1_to_S2(s_network_dir, part1, part2) # Not the fastest, but not too bad.

                    
                    time_constraint_beg = time()

                    cut_s_edges = get_edges_from_S1_to_S2(s_network_dir, part1, part2) # Not the fastest, but not too bad.

                    cut = @build_constraint(sum(y[v_edge, s_edge] for s_edge in cut_s_edges) + sum( x[dst(v_edge), s_other_node] for s_other_node in part1)
                        >=
                        x[src(v_edge), s_node] + x[dst(v_edge), s_node]
                    ) 

                    MOI.submit(model, MOI.UserCut(cb_data), cut)
                    time_constraint += time() - time_constraint_beg


                    #println("Size of cut: $(length(part1))")
                    

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
                nb_cuts = nb_cuts_overall,
                nb_call = nb_call,
                nb_min_cut = nb_min_cut,
                nb_continuity=nb_continuity
        )
    end




    # Get the solution

    
    

    println("Find the solution $(objective_value(model)) in $(time() - time_start)")
    println("Stats branch and cut: $time_seperation_overall s overall, $time_min_cut on mincut, $time_get_values getting values, $time_matrix_flow on flow matrix, $time_constraint on constraint
                $nb_call calls, $nb_min_cut mincuts computed, $nb_cuts_overall cuts found")


    return    ( sol_value= objective_value(model),
                lower_bound = objective_bound(model),
                gap = relative_gap(model),
                node_count = node_count(model),
                time_solving = (time() - time_start),
                time_seperation = time_seperation_overall,
                nb_cuts = nb_cuts_overall, 
                nb_call = nb_call,
                nb_min_cut = nb_min_cut,
                nb_continuity=nb_continuity
    )
end



# I think the conclusion I got is:
# You generate, for cheap, a looooot of cuts.
# Some are useful, lot are useless probably.
# But it still helps the solver. It shows that adding many cuts can be nice in some way...
# I'm surprised by how much it does not slow the solver. But again, it's probably because the useless cuts get removed.
function solve_branch_and_cut_false(instance; solver="cplex", time_solver = 100, alpha_acceptance=0.05)
    
    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    time_start = time()
    if solver == "cplex"
        model = Model(CPLEX.Optimizer)
        set_up_problem_ff_plus(instance, model)
    elseif solver == "gurobi"
        model = Model(Gurobi.Optimizer)
        set_up_problem_ff_plus(instance, model)
        set_optimizer_attribute(model, "Threads", 1)
    else
        println("Solver didnt found. What?")
        return
    end

    set_time_limit_sec(model, time_solver)


    #set_optimizmer_attribute(model, "CPXPARAM_MIP_Strategy_Search", 1)
    #set_optimizer_attribute(model, "CPX_PARAM_REPEATPRESOLVE", 0)
    #
    #set_optimizer_attribute(model, "CPXPARAM_Parallel", -1)

    # setting things for the branch and cut
    n_s = nv(s_network_dir)
    augmented_network = copy(s_network_dir.graph)
    add_vertex!(augmented_network) # n_s + 1: t(\ebar)
    for i in 1:n_s
        if s_network[i][:cap] > 0
            add_edge!(augmented_network, i, n_s+1)
        end
    end



    x = model[:x]
    y = model[:y]

    nb_cuts_overall = 0
    nb_min_cut = 0
    nb_call = 0
    time_seperation_overall = 0
    time_min_cut = 0
    time_get_values = 0
    time_matrix_flow = 0
    time_constraint = 0
    
    matrix_flows = zeros(Float32, n_s+1, n_s+1)

    # Variable stuff
    n_x = nv(v_network) * nv(s_network)
    n_y = ne(v_network) *  ne(s_network_dir)
    x_idx = [MOI.VariableIndex(i) for i in 1:n_x]
    y_idx = [MOI.VariableIndex(i) for i in n_x+1:n_x+n_y]
    x_vals = zeros(Float64, n_x)    
    y_vals = zeros(Float64, n_y)
    backend = JuMP.backend(model)


    function find_cutset_cuts(cb_data)
        
        nb_call += 1
        if nb_call % 100 == 0
            println("At call $nb_call ! Did $nb_min_cut mincuts and found $nb_cuts_overall cuts")
        end


        if alpha_acceptance > 0.99 # Ahah don't worry (this is ugly but it's just for the tests)
            return
        end

        time_start_separation = time()
        nb_cuts_current = 0
        time_beg_getting_values = time()
        
        #=
        x_values = callback_value.(cb_data, x)
        y_values = callback_value.(cb_data, y)
        =#

        MOI.get!(x_vals,
            backend.optimizer.model,
            MOI.CallbackVariablePrimal(cb_data),
            x_idx
        )

        MOI.get!(y_vals,
            backend.optimizer.model,
            MOI.CallbackVariablePrimal(cb_data),
            y_idx
        )

        time_get_values += time()-time_beg_getting_values

        i_v_edge = 0
        for v_edge in edges(v_network)

            time_matrix_flow_beg = time()
            matrix_flows .= 0.0  # resets in place

            i_v_edge += 1
            i_s_edge = 0
            for s_edge in edges(s_network_dir)
                i_s_edge += 1
                matrix_flows[src(s_edge), dst(s_edge)] = y_vals[ ne(v_network) * (i_s_edge-1) + i_v_edge]
            end

            #@inbounds matrix_flows .= clamp.(matrix_flows, 0.0, 1.0)
            @inbounds @simd for i in eachindex(matrix_flows)
                v = matrix_flows[i]
                if !(isfinite(v)) || v < 1e-12
                    matrix_flows[i] = 0.0
                elseif v > 1.0
                    matrix_flows[i] = 1.0
                end
            end

            time_matrix_flow += time() - time_matrix_flow_beg

            #println("matrix flow: $matrix_flows")
            
            for s_node in vertices(s_network_dir)
                val = x_vals[nv(v_network) * (s_node-1) + src(v_edge)]
                
                val < alpha_acceptance * 1.5 && continue # reduces the number of separations
                val > 0.9 && continue # prevent bug in the mincut!

                nb_min_cut += 1

                # Finish flow matrix for that node - remember that no flow between u and t(\ebar) here
                for other_s_node in vertices(s_network_dir)
                    if s_node == other_s_node
                        matrix_flows[other_s_node, n_s+1] = 0.
                    else 
                        matrix_flows[other_s_node, n_s+1] = x_vals[nv(v_network) * (s_node-1) + dst(v_edge)]
                    end
                end
                
                # min cut between s_node and n_s+1
                time_beg_cut = time()
                (part1, part2, flow) = GraphsFlows.mincut(augmented_network, s_node, n_s+1, matrix_flows, DinicAlgorithm())
                time_min_cut += time()-time_beg_cut
                
                if flow < val - alpha_acceptance
                    if n_s+1 ∈ part1
                        println("WOWOWOW NS+1 IS IN PART1?? $part1 and $part2 wthhh")
                        println("I mean, I'm doing a cut between $s_node and $(n_s+1)... Is that algorithm stupid?")
                        break
                    end
                    if s_node ∈ part2
                        println("WOWOWOW $s_node IS IN PART2?? $part1 and $part2 wthhh")
                        println("I mean, I'm doing a cut between $s_node and $(n_s+1)... Is that algorithm stupid?")
                        break
                    end
                    
                    #=
                    println("Damn new cut for $v_edge on $s_node: $part1")
                    for s_node_other ∈ part1
                        if s_node_other != s_node
                            println("       with x_$(s_node_other): $(x_vals[nv(v_network) * (s_node_other-1) + dst(v_edge)])")
                        end
                    end
                    =#

                    #=
                    cut_s_edges = get_edges_from_S1_to_S2(s_network_dir, part1, part2) # Not the fastest, but not too bad.

                    time_constraint_beg = time()

                    terms = Vector{MOI.ScalarAffineTerm{Float64}}()
                    # y part
                    for s_edge in cut_s_edges
                        push!(terms, MOI.ScalarAffineTerm(1.0, y[v_edge, s_edge].index))
                    end
                    # x[dst(v_edge), s_other_node] part
                    for s_other_node in part1
                        if s_other_node != s_node
                            push!(terms, MOI.ScalarAffineTerm(1.0, MOI.VariableIndex(nv(v_network) * (s_other_node-1) + dst(v_edge))))
                        end
                    end
                    push!(terms, MOI.ScalarAffineTerm(-1.0, MOI.VariableIndex(nv(v_network) * (s_node-1) + src(v_edge))))
                
                    
                
                    MOI.submit(backend, MOI.UserCut(cb_data),
                        MOI.ScalarAffineFunction(terms, 0.0),
                        MOI.GreaterThan(0.0)
                    )
                    
                    time_constraint += time() - time_constraint_beg
                    =#
                    
                    time_constraint_beg = time()

                    cut_s_edges = get_edges_from_S1_to_S2(s_network_dir, part1, part2) # Not the fastest, but not too bad.

                    cut = @build_constraint(sum(y[v_edge, s_edge] for s_edge in cut_s_edges) + sum( x[dst(v_edge), s_other_node] for s_other_node in part1)
                        >=
                        x[src(v_edge), s_node] + x[dst(v_edge), s_node]
                    ) 

                    MOI.submit(model, MOI.UserCut(cb_data), cut)
                    time_constraint += time() - time_constraint_beg


                    #println("Size of cut: $(length(part1))")

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
                nb_cuts = nb_cuts_overall,
                nb_call = nb_call
        )
    end




    # Get the solution

    
    

    println("Find the solution $(objective_value(model)) in $(time() - time_start)")
    println("Stats branch and cut: $time_seperation_overall s overall, $time_min_cut on mincut, $time_get_values getting values, $time_matrix_flow on flow matrix, $time_constraint on constraint
                $nb_call calls, $nb_min_cut mincuts computed, $nb_cuts_overall cuts found")


    return    ( sol_value= objective_value(model),
                lower_bound = objective_bound(model),
                gap = relative_gap(model),
                node_count = node_count(model),
                time_solving = (time() - time_start),
                time_seperation = time_seperation_overall,
                nb_cuts = nb_cuts_overall, 
                nb_call = nb_call
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








function solve_simple_cplex(instance; time_solver = 100)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir

    time_start = time()
    model = Model(CPLEX.Optimizer)
    set_up_problem_ff_plus(instance, model)

    set_time_limit_sec(model, time_solver)

    # On
    set_optimizer_attribute(model, "CPXPARAM_Threads", 1)    
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Strategy_Search", 1)

    #set_optimizer_attribute(model, "CPX_PARAM_REPEATPRESOLVE", 0)

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
        time_solving = (time() - time_start),
        nb_call=nb_call
    )
end







# ========== CLASSICAL STUFF


function solve_cutting_plane(instance; alpha_acceptance=0.10)

    time_beginning = time()

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir

    time_start = time()
    model = Model(CPLEX.Optimizer)
    set_up_problem_ff_plus(instance, model)
    relax_integrality(model)
    set_silent(model)


    optimize!(model)
    value_ff = round(objective_value(model); digits=3)

    n_s = nv(s_network_dir)
    augmented_network = deepcopy(s_network_dir.graph)
    add_vertex!(augmented_network) # n_s + 1: t(\ebar)
    for i in 1:n_s
        add_edge!(augmented_network, i, n_s+1)
    end

    nb_cuts_overall = 0

    keep_going = true
    iter = 1
    while keep_going

        optimize!(model)

        status = primal_status(model)
        if status != MOI.FEASIBLE_POINT
            println("It's unfeasible!!")
            return (value_ff, 10e9)
        end
        
        println("Iter $iter, value: $(objective_value(model)), nb cuts: $nb_cuts_overall")

        keep_going = false

        x_values = value.(model[:x])
        y_values = value.(model[:y])

        nb_new_cuts = 0
        
        for v_edge in edges(v_network)

            matrix_flows = zeros(n_s+1, n_s+1)
 
            for s_edge in edges(s_network_dir)
                matrix_flows[src(s_edge), dst(s_edge)] = y_values[v_edge, s_edge]
            end

            for s_node in vertices(s_network_dir)
                #x_values[src(v_edge), s_node] < alpha_acceptance && continue
                
                # Finish flow matrix for that node - remember that no flow between u and t(\ebar) here
                for other_s_node in vertices(s_network_dir)
                    if s_node == other_s_node
                        matrix_flows[other_s_node, n_s+1] = 0.
                    else 
                        matrix_flows[other_s_node, n_s+1] = x_values[dst(v_edge), other_s_node]
                    end
                end

                (part1, part2, flow) = GraphsFlows.mincut(augmented_network, s_node, n_s+1, matrix_flows, DinicAlgorithm())

                if flow < x_values[src(v_edge), s_node] - alpha_acceptance
                    #println("Damn! $part1 to $part2: i got $flow, imma win $(x_values[src(v_edge), s_node] - flow) for $v_edge")
                    cut_s_edges = get_edges_from_S1_to_S2(s_network_dir, part1, part2)
                    @constraint(model, sum(model[:y][v_edge, s_edge] for s_edge in cut_s_edges) + sum( model[:x][dst(v_edge), s_other_node] for s_other_node in part1)
                        >=
                        model[:x][src(v_edge), s_node] + model[:x][dst(v_edge), s_node]
                    )

                    nb_cuts_overall += 1
                    keep_going = true
                    nb_new_cuts += 1
                end
            end


        end

        println("I generated $nb_new_cuts new cuts!")
        iter += 1

    end



    value_cutting_plane = round(objective_value(model); digits=3)
    println("\n\n\n FINISHED, obtained $value_cutting_plane in $(time() - time_beginning)s")

    return (value_ff, value_cutting_plane)
end














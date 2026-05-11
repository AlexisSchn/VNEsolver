
using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf
using Statistics


#general
includet("../../utils/import_utils.jl")

# utils colge
includet("utils/utils-subgraphdecompo.jl")
includet("utils/partition-graph.jl")
includet("utils/checkers.jl")

# pricers
includet("pricers/pricer-exact.jl")
includet("pricers/pricer-subsn.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")
includet("end-heuristic/diving-heuristic.jl")
includet("end-heuristic/local-search-exact.jl")



function solve_subgraph_decompo_one_hour(instance)


    # Budget: 3600 seconds
    time_init = 250
    time_limit_sub_pricers = 1800
    time_limit_for_full_pricers = 3300
    time_cg_heuristic = 300
    #time_local_search = 300


    println("Starting...")
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network


    # ======= SETTING UP THE DECOMPOSITION ======= #
    nb_virtual_subgraph = floor(Int, nv(v_network.graph)/10)
    v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.2)

    vn_decompo = set_up_decompo(instance, v_node_partitionning)
    
    println("Decomposition set: ")
        println("For $v_network, there is $(length(vn_decompo.subgraphs)) subgraphs:")

    for subgraph in vn_decompo.subgraphs
        println("       $(subgraph.graph[][:name]) with $(nv(subgraph.graph)) nodes")
    end
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")

    
    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    print("Master problem set... ")



    #------------ GENERATION DE COLONNES
    nb_columns = 0
    nb_iter = 0
    nb_usable = 0

    time_master = 0
    time_subproblems = 0

    cg_value = 99999
    lower_bound =  0   # Computing base linear relaxation ?



    print("\n\n==================== Starting CG ====================\n")



    # ====== STEP 1: init
    println("------- Part 1: Initialization")
    print("Paving the substrate network... ")

    time_beginning_init = time()
    nb_substrate_subgraph = 20
    nb_nodes_subgraph = 25
    sn_decompo_clusters = get_sn_decompo(s_network, nb_substrate_subgraph, nb_nodes_subgraph)

    time_limit_per_pb = time_init * 2. / (nb_substrate_subgraph*length(vn_decompo.subgraphs))
    for vn_subgraph in vn_decompo.subgraphs
        pricers = set_up_pricer_sn_decompo(instance, vn_subgraph, sn_decompo_clusters, "normal")
        nb_col_cur = 0
        for sub_pricer in pricers
            column, true_cost, reduced_cost = solve_pricers_sn_decompo(sub_pricer, time_limit=time_limit_per_pb)
            if column !== nothing && true_cost < 9999999
                add_column(master_problem, instance, sub_pricer.vn_subgraph, column, true_cost)
                nb_columns += 1
                nb_col_cur+=1
            end
        end
        print("done for $(vn_subgraph.graph[][:name]), $nb_col_cur columns found; ")
    end
    println("Initialization complete, $nb_columns columns found, in $(time() - time_beginning_init)\n")
    # First iteration of master problem
    optimize!(model)
    time_master +=  solve_time(model)
    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("Infeasible or unfinished: $status")
        return
    end

    cg_value = objective_value(model)
    time_overall = time() - time_beginning

    @printf("Iter %2d  CG bound: %10.3f  lower bound: %10.3f  %5d column  time: %5.2fs  \n\n",
        nb_iter, cg_value, lower_bound, nb_columns, time_overall
    )




    # ====== STEP 2: smaller pricers - paving the network
    println("------- Part 2: Reduced pricers")

    time_beginning_sub_pricers = time()

    nb_substrate_subgraph = floor(Int, nv(s_network) / 15)  
    nb_nodes_subgraph = 30
    sn_decompo_clusters = get_sn_decompo(s_network, nb_substrate_subgraph, nb_nodes_subgraph)
    println("We have $nb_substrate_subgraph sub-substrate, with at least $nb_nodes_subgraph capacited nodes")

    pricers_sn_decompo = OrderedDict()
    sub_pricers_last_values = OrderedDict()
    nb_desactivated_pricers = 0
    desactivated_pricers = []     
    nb_sub_pricers = nb_substrate_subgraph * length(vn_decompo.subgraphs)
    for vn_subgraph in vn_decompo.subgraphs
        pricers_sn_decompo[vn_subgraph] = set_up_pricer_sn_decompo(instance, vn_subgraph, sn_decompo_clusters, "normal")
        for pricer in pricers_sn_decompo[vn_subgraph]
            sub_pricers_last_values[pricer] = -99999.
        end
    end

    keep_on = true
    reason = "I don't know"
    while keep_on
        nb_iter += 1

        # ---- pricers part

        dual_costs = get_duals(instance, vn_decompo, master_problem)
        
        for key in keys(sub_pricers_last_values)
            sub_pricers_last_values[key] = sub_pricers_last_values[key]*1.1
        end
        
        sorted_subpb = sort(collect(sub_pricers_last_values), by=x->x[2])
        average_obj = 0
        nb_pricer_to_do = min(5, length(keys(sub_pricers_last_values)))
        for couple in sorted_subpb[1:nb_pricer_to_do]


            pricer_sub_sn = couple[1]

            time_limit_subpb = minimum([10., time_limit_sub_pricers - (time()-time_beginning_sub_pricers)])
            if time_limit_subpb < 0.01
                break
            end

            update_pricer_sn_decompo(vn_decompo, pricer_sub_sn, dual_costs)
            column, true_cost, reduced_cost = solve_pricers_sn_decompo(pricer_sub_sn, time_limit=time_limit_subpb)

            if column !== nothing && reduced_cost < -0.001
                add_column(master_problem, instance, pricer_sub_sn.vn_subgraph, column, true_cost)
                nb_columns += 1
            end

            if reduced_cost>-0.001
                if pricer_sub_sn ∉ desactivated_pricers
                    nb_desactivated_pricers+=1
                    push!(desactivated_pricers, pricer_sub_sn)
                end
            end
            sub_pricers_last_values[pricer_sub_sn] = reduced_cost

            average_obj += (reduced_cost/nb_pricer_to_do)
            time_subproblems += solve_time(pricer_sub_sn.model) 
        end


        # ---- master problem part

        optimize!(model)
        time_master +=  solve_time(model)

        cg_value = objective_value(model)

        time_overall = time()-time_beginning

        @printf("Iter %2d  CG bound: %10.3f  lower bound: %10.3f  %5d column  time: %5.2fs  average reduced cost: %10.3f \n",
                    nb_iter, cg_value, lower_bound, nb_columns, time_overall, average_obj)
 


        
        # ----- useful things

        time_sub_pricers = time() - time_beginning_sub_pricers
        if time_sub_pricers < time_limit_sub_pricers
            keep_on = true
            if nb_desactivated_pricers >= nb_sub_pricers   
                keep_on = false
                reason="changing to full pricers to get better columns"
            end
            if nb_columns>1300 #nb_columns>250*length(vn_decompo.subgraphs) 
                keep_on=false
                reason="too many columns generated already..."
            end
        else
            keep_on = false
            reason="time limit"
        end

    end
    println("\n Step 2 finished, reason: $reason.")
    


    # ====== STEP 3: full pricers
    println("\n------- Solving method: Exact pricers")
    time_for_full_pricers = time_limit_for_full_pricers - (time() -time_beginning) 
    time_beginning_full_pricers = time()

    pricers_full = Dict()
    for subgraph in vn_decompo.subgraphs
        #pricers_full[subgraph] = set_up_pricer_cons(instance, subgraph)
        pricers_full[subgraph] = set_up_pricer(instance, subgraph)
    end
    
    keep_on = true
    reason = "I don't know"
    while keep_on
        nb_iter += 1

        # ---- Pricers things
        dual_costs = get_duals(instance, vn_decompo, master_problem)

        has_found_new_column = false
        has_finished = true
        sum_pricers_values = 0

        for vn_subgraph in vn_decompo.subgraphs
            pricer = pricers_full[vn_subgraph]

            time_limit_pricer = time_for_full_pricers - (time()-time_beginning_full_pricers)
            if time_limit_pricer < 0.1
                has_finished = false
                break
            end

            sub_mapping, true_cost, reduced_cost = update_solve_pricer(instance, vn_decompo, pricer, dual_costs; time_limit = time_limit_pricer)

            
            if (!isnothing(sub_mapping)) && reduced_cost < -0.0001
                has_found_new_column = true
                add_column(master_problem, instance, vn_subgraph, sub_mapping, true_cost)
                nb_columns += 1
            end
            
            if isnothing(sub_mapping)
                println("Pricer with no solution found, stopping the CG...")
                reason="pricer-unfeasible"
                has_found_new_column = false
                break
            end

            sum_pricers_values += reduced_cost
            time_subproblems += solve_time(pricer.model) 

            time_limit_pricer = time_for_full_pricers - (time()-time_beginning_full_pricers)
            if time_limit_pricer < 0.1
                has_finished = false
                break
            end
        end

        if has_finished
            current_lower_bound = cg_value + sum_pricers_values
            if current_lower_bound > lower_bound && nb_iter > 3
                lower_bound = current_lower_bound
            end
        end


        # ----- Master problem stuff

        optimize!(model)
        time_master +=  solve_time(model)

        cg_value = objective_value(model)

        time_overall = time()-time_beginning
        average_obj = sum_pricers_values/length(vn_decompo.subgraphs)

        @printf("Iter %2d  CG bound: %10.3f  lower bound: %10.3f  %5d column  time: %5.2fs  average reduced cost: %10.3f \n",
            nb_iter, cg_value, lower_bound, nb_columns, time_overall, average_obj)


        if (time() - time_beginning_full_pricers) < time_for_full_pricers 
            keep_on = true
            if !has_found_new_column
                keep_on = false
                reason="no improving columns"
            end
        else
            keep_on = false
            reason="time limit"
        end


    end

        


    print("\n==================== CG finished ====================\nReason: $reason \n")
    println("Time in MP: $(round(time_master; digits=3)) , time in SP: $(round(time_subproblems; digits=3)), time overall: $(round(time_overall; digits=3))")
    println("$nb_iter iters, final value: $(round(cg_value; digits=3))")
    println("====================================================\n")




    # ======= END HEURISTICS ======= #

    # ---- Price n Branch heuristic
    time_cg_heuristic = 300
    value_cg_heuristic, cg_heuristic_solution = basic_heuristic(instance, vn_decompo, master_problem, time_cg_heuristic)


    #= ---- Large Neighbor Search
    #local_search(instance, vn_decompo, heur_sol)
    if cg_heuristic_solution != nothing
        value_local_search, lsn_placement_sol = local_search_changin(instance, cg_heuristic_solution, time_local_search)
    else
        value_local_search = 9999
    end
    =#

    result = Dict()
    result["algorithm"] = "one-hour"
    result["time"] = time() - time_beginning_init
    result["cg_value"] = cg_value
    result["lower_bound"] = lower_bound
    result["nb_iter"] = nb_iter
    result["nb_col"] = nb_columns
    result["value_cg_heuristic"] = value_cg_heuristic
    #result["value_local_search"] = value_local_search


    return result
end



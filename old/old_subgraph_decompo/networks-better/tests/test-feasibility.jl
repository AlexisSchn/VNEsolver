
using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf
using Statistics

#general
includet("../../../utils/import_utils.jl")

# utils colge
includet("../utils/utils-subgraphdecompo.jl")
includet("../utils/partition-vn.jl")
includet("../utils/checkers.jl")

# pricers
includet("../pricers/pricer-exact.jl")
includet("../pricers/pricer-subsn.jl")

# end heuristics
includet("../end-heuristic/basic-ilp.jl")



function solve_subgraph_decompo(instance; time_max = 3600, v_node_partitionning = [], nb_part = 2, type_pricer="normal")

    println("Starting...")
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network


    # ======= SETTING UP THE DECOMPOSITION ======= #
    if v_node_partitionning == []
        if nb_part<0
            nb_part = floor(Int, nv(v_network.graph)/9)+1
        end
        v_node_partitionning = partition_vn(instance, nb_part)
    end

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
    time_overall = time()-time_beginning
    part_colge = 0.90
    time_colge = part_colge * time_max 
    time_end_solving = (1-part_colge) * time_max 

    cg_value = 99999
    lower_bound =  0   # Computing base linear relaxation ?

    print("\n\n==================== Starting CG ====================\n")


    # ====== STEP 2: smaller pricers - paving the network
    println("------- Part 2: Reduced pricers")

    nb_parts = 3
    nb_nodes_per_snsubgraph_pricers = 25
    println("$nb_parts substrate subgraphs to do, with at least $nb_nodes_per_snsubgraph_pricers capacitated nodes...")
    sn_decompo_clusters = get_sn_decompo_kahip(s_network, nb_parts, nb_nodes_per_snsubgraph_pricers)
    pricers_sn_decompo = OrderedDict()
    for vn_subgraph in vn_decompo.subgraphs
        pricers_sn_decompo[vn_subgraph] = set_up_pricer_sn_decompo(instance, vn_subgraph, sn_decompo_clusters, type_pricer)
    end

    keep_on = true
    reason = "I don't know"
    while keep_on
        nb_iter += 1

        # ---- master problem
        optimize!(model)
        time_master +=  solve_time(model)
        status = termination_status(model)
        if status != MOI.OPTIMAL
            println("Infeasible or unfinished: $status")
            return
        end

        cg_value = objective_value(model)
        dual_costs = get_duals(instance, vn_decompo, master_problem)

        # ---- pricers
        time_beginning_pricer = time()


        average_obj = 0

        has_found_new_column = false
        sum_pricers_values = 0
        for vn_subgraph in vn_decompo.subgraphs
            for pricer in pricers_sn_decompo[vn_subgraph]
                

                update_pricer_sn_decompo(vn_decompo, pricer, dual_costs)
                column, true_cost, reduced_cost = solve_pricers_sn_decompo(pricer, time_limit=500)

    
                if column !== nothing && reduced_cost < -0.01
                    has_found_new_column = true
                    add_column(master_problem, instance, vn_subgraph, column, true_cost)
                    nb_columns += 1
                end
    
                sum_pricers_values += reduced_cost
            end
        end
        average_obj = sum_pricers_values / (nb_parts * length(vn_decompo.subgraphs))
        

        @printf("Iter %2d  CG bound: %10.3f  lower bound: %10.3f  %5d column  time: %5.2fs  average reduced cost: %10.3f \n",
                    nb_iter, cg_value, lower_bound, nb_columns, time_overall, average_obj)
    

        time_subproblems += time() - time_beginning_pricer
        time_overall = time()-time_beginning
        keep_on = false
        if time_overall < time_colge
            if has_found_new_column
                keep_on = true
            end
        else
            keep_on = false
            reason="time limit"
        end

    end
    println("\n Step 2 finished, reason: $reason.")


        

    # ====== STEP 3: full pricers
    println("\n------- Solving method: Exact pricers")
    reason = "I don't know"

    #=
    pricers_full = Dict()
    for subgraph in vn_decompo.subgraphs
        if type_pricer == "ghost"
            pricers_full[subgraph] = set_up_pricer_ghost(instance, subgraph)
        elseif type_pricer == "constraint"
            pricers_full[subgraph] = set_up_pricer_cons(instance, subgraph)
        else
            pricers_full[subgraph] = set_up_pricer(instance, subgraph)
        end
    end
    keep_on = true
    reason = "I don't know"
    while keep_on
        nb_iter += 1

        optimize!(model)
        time_master +=  solve_time(model)

        status = termination_status(model)
        if status != MOI.OPTIMAL
            println("Infeasible or unfinished: $status")
            return
        end
        cg_value = objective_value(model)

        dual_costs = get_duals(instance, vn_decompo, master_problem)


        has_found_new_column = true
        time_beginning_pricer = time()

        # do all subpb, get lg bound ?
        # also need to ensure that time is respected !
        has_found_new_column = false
        sum_pricers_values = 0
        for vn_subgraph in vn_decompo.subgraphs
            pricer = pricers_full[vn_subgraph]

            time_limit_subpb = time_colge - (time()-time_beginning)
            if time_limit_subpb < 0.1
                break
            end

            sub_mapping, true_cost, reduced_cost = update_solve_pricer(instance, vn_decompo, pricer, dual_costs; time_limit = time_limit_subpb)

            if reduced_cost < -0.0001
                has_found_new_column = true
                add_column(master_problem, instance, vn_subgraph, sub_mapping, true_cost)
                nb_columns += 1
            end

            sum_pricers_values += reduced_cost
        end

        time_limit_subpb = time_colge - (time()-time_beginning)
        if time_limit_subpb > 0.1 # computing the lagrangian bound, only if I've solved all subpb
            current_lower_bound = cg_value + sum_pricers_values
            if current_lower_bound > lower_bound && nb_iter > 3
                lower_bound = current_lower_bound
            end
        end

        time_subproblems += time() - time_beginning_pricer

        average_obj = sum_pricers_values/length(vn_decompo.subgraphs)
        @printf("Iter %2d  CG bound: %10.3f  lower bound: %10.3f  %5d column  time: %5.2fs  average reduced cost: %10.3f \n",
            nb_iter, cg_value, lower_bound, nb_columns, time_overall, average_obj)



        time_overall = time()-time_beginning
        if time_overall < time_max * 0.90
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
    =#

    print("\n==================== CG finished ====================\nReason: $reason \n")
    println("Time in MP: $(round(time_master; digits=3)) , time in SP: $(round(time_subproblems; digits=3)), time overall: $(round(time_overall; digits=3))")
    println("$nb_iter iters, final value: $(round(cg_value; digits=3))")
    println("====================================================\n")

    check_all_columns(vn_decompo, master_problem)
    
    nb_valid_cols, values = check_all_columns_2(vn_decompo, master_problem, instance)
    if nb_valid_cols > 0
        mean_vals = mean(values)
        min_vals = minimum(values)
        median_vals = median(values)
    else
        mean_vals = 0
        min_vals = 0
        median_vals = 0
    end
    # ======= END HEURISTIC STUFF ======= #

    time_end_solving = 6000
    val, heur_sol = basic_heuristic(instance, vn_decompo, master_problem, time_end_solving)

    result = Dict()
    result["time"] = round(time() - time_beginning, digits=2)
    result["CG_val"] = round(cg_value, digits=2)
    result["LG_val"] = round(lower_bound, digits=2)
    result["nb_iter"] = nb_iter
    result["nb_col"] = nb_columns
    result["heuristic_res"] = round(val, digits=2)
    result["nb_valid_cols"] = nb_valid_cols
    result["mean_vals"] = round(mean_vals, digits=2)
    result["min_vals"] = round(min_vals, digits=2)
    result["median_vals"] = round(median_vals, digits=2)

    return result
end



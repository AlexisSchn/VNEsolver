# todo for sauron
# => pave the network with cheapest heuristics
# => remove all the useless stuff, make the colge as clean as possible for ifip
# => better paving with cplex : get the right sizes that work better, get several solution at the same time


using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf

using DataFrames
using CSV

#general
includet("../../../utils/import_utils.jl")
#includet("../../utils/visu.jl")

# utils colge
includet("../utils/utils-subgraphdecompo.jl")
includet("../utils/partition-vn.jl")
includet("../utils/base_relaxation.jl")

# pricers
includet("../pricers/pricer-full.jl")
includet("../pricers/sn-decompo.jl")

# end heuristics
includet("../end-heuristic/basic-ilp.jl")



function solve_subgraph_decompo(instance; time_max = 100, v_node_partitionning = [], nb_part = -1)

    df_stats = DataFrame(Time=Float32[], CG_val=Float32[], LB_val=Float32[], Avg_red=Float32[])




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

    time_master = 0
    time_subproblems = 0
    time_overall = time()-time_beginning
    part_colge = 0.90
    time_colge = part_colge * time_max 
    time_end_solving = (1-part_colge) * time_max 

    cg_value = 99999
    lower_bound =  0   # Computing base linear relaxation ?

    print("\n\n==================== Starting CG ====================\n")


    # ====== STEP 0: computing a lower bound, through flow formulation.
    println("-------Part 0: lower bound through flow formulation")
    lower_bound = get_base_relaxation(instance)
    println("Bound obtained: $lower_bound\n\n")

    #=
    # ====== STEP 1: paving with ILP
    println("------- Part 1: Initialization")
    print("Paving the substrate network... ")
    dual_costs = get_empty_duals(instance, vn_decompo)
    nb_max_node_subgraphs = maximum([nv(subgraph.graph) for subgraph in vn_decompo.subgraphs])
    nb_parts_init = floor(Int, nv(s_network) / (nb_max_node_subgraphs))+1
    nb_nodes_per_snsubgraph_init = max(15, 1.5 * nb_max_node_subgraphs)
    sn_decompo_clusters_init = get_sn_decompo(s_network, nb_parts_init, nb_nodes_per_snsubgraph_init)
    init_problems = Dict()
    for vn_subgraph in vn_decompo.subgraphs
        init_problems[vn_subgraph] = set_up_pricer_sn_decompo(instance, vn_subgraph, sn_decompo_clusters_init)
        nb_col_cur = 0
        for problem in init_problems[vn_subgraph]
            update_pricer_sn_decompo(vn_decompo, problem, dual_costs) #ugh
            column, true_cost, reduced_cost = solve_pricers_sn_decompo(problem, time_limit=10)
            if column !== nothing && true_cost < 9999999
                add_column(master_problem, instance, problem.vn_subgraph, column, true_cost)
                nb_columns += 1
                nb_col_cur+=1
            end
        end
        print("done for $(vn_subgraph.graph[][:name]), $nb_col_cur columns found; ")
    end
    println("Initialization complete, $nb_columns columns found\n\n")
  =#






    # ====== STEP 2: still some paving: smaller SN pricers 
    println("------- Part 2: Reduced pricers")

    nb_max_node_subgraphs = maximum([nv(subgraph.graph) for subgraph in vn_decompo.subgraphs])

    nb_parts = floor(Int, nv(s_network) / (nb_max_node_subgraphs)/2) + 1
    nb_nodes_per_snsubgraph_pricers = max(25, 2.5 * nb_max_node_subgraphs)
    println("$nb_parts substrate subgraphs to do, with at least $nb_nodes_per_snsubgraph_pricers capacitated nodes...")
    sn_decompo_clusters = get_sn_decompo(s_network, nb_parts, nb_nodes_per_snsubgraph_pricers)
    pricers_sn_decompo = OrderedDict()
    sn_decompo_last_iter_sol = OrderedDict()
    nb_desactivated_pricers = 0
    desactivated_pricers = []     
    nb_pricers = nb_parts * length(vn_decompo.subgraphs)
    for vn_subgraph in vn_decompo.subgraphs
        pricers_sn_decompo[vn_subgraph] = set_up_pricer_sn_decompo(instance, vn_subgraph, sn_decompo_clusters)
        for pricer in pricers_sn_decompo[vn_subgraph]
            sn_decompo_last_iter_sol[pricer] = -99999.
        end
    end

    keep_on = false
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
        
        for key in keys(sn_decompo_last_iter_sol)
            sn_decompo_last_iter_sol[key] = sn_decompo_last_iter_sol[key]*1.1
        end
        

        sorted_subpb = sort(collect(sn_decompo_last_iter_sol), by=x->x[2])
        average_obj = 0
        nb_pricer_to_do = min(5, length(keys(sn_decompo_last_iter_sol)))
        for couple in sorted_subpb[1:nb_pricer_to_do]


            pricer_sub_sn = couple[1]


            time_limit_subpb = time_colge - (time()-time_beginning)
            if time_limit_subpb < 0.01
                break
            end

            update_pricer_sn_decompo(vn_decompo, pricer_sub_sn, dual_costs)
            column, true_cost, reduced_cost = solve_pricers_sn_decompo(pricer_sub_sn, time_limit=time_limit_subpb)

            if column !== nothing && reduced_cost < -0.001
                add_column(master_problem, instance, pricer_sub_sn.vn_subgraph, column, true_cost)
                nb_columns += 1
            end

            if reduced_cost>-4.
                if pricer_sub_sn ∉ desactivated_pricers
                    nb_desactivated_pricers+=1
                    push!(desactivated_pricers, pricer_sub_sn)
                end
            end
            sn_decompo_last_iter_sol[pricer_sub_sn] = reduced_cost

            average_obj += (reduced_cost/nb_pricer_to_do)
        end

        @printf("Iter %2d  CG bound: %10.3f  lower bound: %10.3f  %5d column  time: %5.2fs  average reduced cost: %10.3f \n",
                    nb_iter, cg_value, lower_bound, nb_columns, time_overall, average_obj)
 

        time_subproblems += time() - time_beginning_pricer
        time_overall = time()-time_beginning
        if time_overall < time_colge
            keep_on = true
            if nb_desactivated_pricers >= nb_pricers   
                keep_on = false
                reason="changing to full solving to get better columns"
            end 
            if nb_columns>150*length(vn_decompo.subgraphs)
                keep_on=false
                reason="too many columns generated already..."
            end
        else
            keep_on = false
            reason="time limit"
        end


        push!(df_stats, (time_overall, cg_value, lower_bound, average_obj)) 



    end
    println("\n Step 2 finished, reason: $reason. By the way: there was $nb_pricers")




    CSV.write("decompo_stats.csv", df_stats)


    # ====== STEP 3: full pricers
    println("\n\n\n------- Part 3: Exact pricers")

    pricers_full = Dict()
    for subgraph in vn_decompo.subgraphs
        pricers_full[subgraph] = set_up_pricer(instance, subgraph)
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

            column, true_cost, reduced_cost = update_solve_pricer(instance, vn_decompo, pricer, dual_costs; time_limit = time_limit_subpb)

            if reduced_cost < -0.0001
                has_found_new_column = true 
                add_column(master_problem, instance, vn_subgraph, column, true_cost)
                nb_columns += 1
            end

            sum_pricers_values += reduced_cost
        end

        time_limit_subpb = time_colge - (time()-time_beginning)
        if time_limit_subpb > 0.1 # computing the lagrangian bound, only if I've solved all subpb
            current_lower_bound = cg_value + sum_pricers_values
            if current_lower_bound > lower_bound
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

        push!(df_stats, (time_overall, cg_value, lower_bound, average_obj)) 


    end

        


    print("\n==================== CG finished ====================\nReason: $reason \n")
    println("Time in MP: $(round(time_master; digits=3)) , time in SP: $(round(time_subproblems; digits=3)), time overall: $(round(time_overall; digits=3))")
    println("$nb_iter iters, final value: $(round(cg_value; digits=3))")
    println("====================================================\n")


    CSV.write("decompo_stats.csv", df_stats)


    # ======= END HEURISTIC STUFF ======= #

    solution_heuristic = basic_heuristic(instance, vn_decompo, master_problem, time_end_solving)

    return solution_heuristic
    
end

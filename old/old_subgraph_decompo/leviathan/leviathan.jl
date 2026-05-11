# LEVIATHAN
# Trying to be better in the long run.
# Trying to be better in the SN decomposition.
# Tried automatic vn_decompo but it didnt work very well.
# Also, I tried to be a bit safer, this is "production ready".
# Also, it should run lagrangian bound correctly

using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections

#general
includet("../../utils/import_utils.jl")
#includet("../../utils/visu.jl")

# utils colge
includet("utils/utils-subgraphdecompo.jl")
includet("utils/partition-vn.jl")

# pricers
includet("pricers/pricer-full.jl")
includet("pricers/sn-decompo.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")



function solve_subgraph_decompo(instance; time_max = 60, v_node_partitionning = [], nb_part = -1)

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
        println("For $v_network, there is: 
            $(length(vn_decompo.subgraphs)) subgraphs, 
            $(length(vn_decompo.v_nodes_master)) nodes in no subgraph,
            $(length(vn_decompo.v_edges_master)) edges in no subgraph")

    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    print("Master problem set... ")





    # ======= SETTING UP THE PRICERS ======= #
    dual_costs = get_empty_duals(instance, vn_decompo)

    nb_max_node_subgraphs = maximum([nv(subgraph.graph) for subgraph in vn_decompo.subgraphs])
    #----- sn partitionning smaller pricers
    nb_part_smaller = floor(Int, nv(s_network) / (nb_max_node_subgraphs))
    smaller_sn_decompo_clusters = get_sn_decompo(s_network, nb_part_smaller, nb_max_node_subgraphs*1.5)
    pricers_smaller_sn_decompo = OrderedDict()
    smaller_sn_decompo_last_iter_sol = OrderedDict()
    smaller_sn_decompo_nb_desactivated = 0
    smaller_sn_decompo_overall_subpb = nb_part_smaller * length(vn_decompo.subgraphs)
    some_i = 1
    for vn_subgraph in vn_decompo.subgraphs
        pricers_smaller_sn_decompo[vn_subgraph] = set_up_pricer_sn_decompo(instance, vn_subgraph, smaller_sn_decompo_clusters)
        for pricer in pricers_smaller_sn_decompo[vn_subgraph]
            smaller_sn_decompo_last_iter_sol[pricer] = -99999. -some_i
            some_i += 1
        end
    end

    #----- sn partittioning bigger pricers
    nb_part_bigger = floor(Int, nv(s_network) / (nb_max_node_subgraphs*2.75))
    bigger_sn_decompo_clusters = get_sn_decompo(s_network, nb_part_bigger, nb_max_node_subgraphs*3.5)
    pricers_bigger_sn_decompo = OrderedDict()
    bigger_sn_decompo_last_iter_sol = OrderedDict()
    bigger_sn_decompo_nb_desactivated = 0
    bigger_sn_decompo_overall_subpb = nb_part_bigger * length(vn_decompo.subgraphs)
    some_i = 1
    for vn_subgraph in vn_decompo.subgraphs
        pricers_bigger_sn_decompo[vn_subgraph] = set_up_pricer_sn_decompo(instance, vn_subgraph, bigger_sn_decompo_clusters)
        for pricer in pricers_bigger_sn_decompo[vn_subgraph]
            bigger_sn_decompo_last_iter_sol[pricer] = -99999. -some_i
            some_i += 1
        end
    end


    #------ classical pricers
    pricers_full = Dict()
    for subgraph in vn_decompo.subgraphs
        pricers_full[subgraph] = set_up_pricer(instance, subgraph)
    end



    #------------ GENERATION DE COLONNES
    nb_columns = 0
    nb_iter = 0

    time_master = 0
    time_subproblems = 0
    time_overall = time()-time_beginning
    part_colge = 0.90
    time_colge = part_colge * time_max 
    time_end_solving = (1-part_colge) * time_max 

    keep_on = true
    algo_to_use = 1 # 1 is SN partitionning, 2 is reduced pricers, 3 is full pricers
    reason = "I don't know"
    CG_bound = 99999
    LG_bound = 0


    #----- Generating first columns => todo with heuristic in next round
    for vn_subgraph in vn_decompo.subgraphs
        for pricer_sub_sn in pricers_smaller_sn_decompo[vn_subgraph]
            update_pricer_sn_decompo(vn_decompo, pricer_sub_sn, dual_costs)
            column, obj_value = solve_pricers_sn_decompo(pricer_sub_sn)
            print(".")
            if obj_value < 99999
                add_column(master_problem, instance, vn_subgraph, column)
                nb_columns += 1
            end
        end
    end

    print("First set of $nb_columns columns generated, ")



    print("\n\n==================== Starting CG ====================\n")



    while keep_on
        nb_iter += 1
        print("Iter $nb_iter, ")

        optimize!(model)
        time_master +=  solve_time(model)

        status = termination_status(model)
        if status != MOI.OPTIMAL
            println("Infeasible or unfinished: $status")
            return
        end
        CG_bound = objective_value(model)
        print("CG value : $(round(CG_bound;   digits=3)), ")

        dual_costs = get_duals(instance, vn_decompo, master_problem)



        # ============ PRICER STUFF ============
        
        has_found_new_column = true
        time_beginning_pricer = time()

        # ------- SMALLER SN PARTITION ---------
        if algo_to_use == 1  
            print("pricer=smallerpaving, ")

            sorted_subpb = sort(collect(smaller_sn_decompo_last_iter_sol), by=x->x[2])
            overall_obj = 0


            for couple in sorted_subpb[1:5]
                pricer_sub_sn = couple[1]
                time_limit_subpb = time_colge - (time()-time_beginning)
                if time_limit_subpb < 0.01
                    break
                end

                update_pricer_sn_decompo(vn_decompo, pricer_sub_sn, dual_costs)
                column, obj_value = solve_pricers_sn_decompo(pricer_sub_sn, time_limit=time_limit_subpb)

                if column !== nothing
                    add_column(master_problem, instance, pricer_sub_sn.vn_subgraph, column)
                    nb_columns += 1
                end

                if obj_value>-5.
                    smaller_sn_decompo_nb_desactivated+=1
                end
                smaller_sn_decompo_last_iter_sol[pricer_sub_sn] = obj_value

                overall_obj += obj_value
            end
            print("mean sol: $(round(overall_obj/5;digits=3)), ")


            if smaller_sn_decompo_nb_desactivated >= smaller_sn_decompo_overall_subpb   
                print("too many subpb desactivated - changing pricers")
                algo_to_use=2
            end 


        elseif algo_to_use == 2  
            print("pricer=biggerpaving, ")

            sorted_subpb = sort(collect(bigger_sn_decompo_last_iter_sol), by=x->x[2])
            overall_obj = 0
            for couple in sorted_subpb[1:5]
                pricer_sub_sn = couple[1]

                time_limit_subpb = time_colge - (time()-time_beginning)
                if time_limit_subpb < 0.01
                    break
                end

                update_pricer_sn_decompo(vn_decompo, pricer_sub_sn, dual_costs)
                column, obj_value = solve_pricers_sn_decompo(pricer_sub_sn, time_limit=time_limit_subpb)

                if column !== nothing
                    add_column(master_problem, instance, pricer_sub_sn.vn_subgraph, column)
                    nb_columns += 1
                end

                if obj_value>-5.
                    bigger_sn_decompo_nb_desactivated+=1
                end
                bigger_sn_decompo_last_iter_sol[pricer_sub_sn] = obj_value

                overall_obj += obj_value
            end
            print("mean sol: $(round(overall_obj/5;digits=3)), ")


            if bigger_sn_decompo_nb_desactivated >= bigger_sn_decompo_overall_subpb   
                print("too many subpb desactivated - changing pricers")
                algo_to_use=3
            end 


                
        # --------- FULL PRICING -----------
        elseif algo_to_use == 3 # FULL PRICING
            print("exact, ")

            # do all subpb, get lg bound ?
            # also need to ensure that time is respected !
            has_found_new_column = false
            sum_pricers_values = 0
            for vn_subgraph in vn_decompo.subgraphs
                pricer = pricers_full[vn_subgraph]

                time_limit_subpb = time_colge - (time()-time_beginning)
                if time_limit_subpb < 0.01
                    break
                end

                column, obj_value = update_solve_pricer(instance, vn_decompo, pricer, dual_costs; time_limit = time_limit_subpb)

                if obj_value < -0.0001
                    has_found_new_column = true 
                    add_column(master_problem, instance, vn_subgraph, column)
                    nb_columns += 1
                    push!(vn_subgraph.columns, column)
                end

                sum_pricers_values += obj_value
            end

            # need to compute the LG bound here !
            LG_bound = CG_bound + sum_pricers_values
            print("LG bound: $(round(LG_bound; digits=3)), ")

        else
            println("pricer: not cooked yet! ")
            return
        end
        time_subproblems += time() - time_beginning_pricer

        print("with $nb_columns cols")



        time_overall = time()-time_beginning
        if time_overall < time_max * 0.90
            keep_on = true
            if algo_to_use == 2
                if !has_found_new_column
                    keep_on = false
                    reason="no improving columns"
                end
            end

        else
            keep_on = false
            reason="time limit"
        end

        print("in $(round(time_overall;digits=2))s")



        println("")
    end


        


    print("\n==================== CG finished ====================\nReason: $reason \n")
    println("Time in MP: $(round(time_master; digits=3)) , time in SP: $(round(time_subproblems; digits=3)), time overall: $(round(time_overall; digits=3))")
    println("$nb_iter iters, final value: $(round(CG_bound; digits=3))")
    println("====================================================\n")



    # ======= END HEURISTIC STUFF ======= #

    solution_heuristic = basic_heuristic(instance, vn_decompo, time_end_solving)

    return solution_heuristic
    
end




function solve_subgraph_decompo_tests(instance; time_max = 60, v_node_partitionning = [], nb_part = -1)

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
        println("For $v_network, there is: 
            $(length(vn_decompo.subgraphs)) subgraphs, 
            $(length(vn_decompo.v_nodes_master)) nodes in no subgraph,
            $(length(vn_decompo.v_edges_master)) edges in no subgraph")

    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    print("Master problem set... ")





    # ======= SETTING UP THE PRICERS ======= #
    dual_costs = get_empty_duals(instance, vn_decompo)

    nb_max_node_subgraphs = maximum([nv(subgraph.graph) for subgraph in vn_decompo.subgraphs])
    #----- sn partitionning smaller pricers
    nb_part_smaller = floor(Int, nv(s_network) / (nb_max_node_subgraphs))
    smaller_sn_decompo_clusters = get_sn_decompo(s_network, nb_part_smaller, nb_max_node_subgraphs*1.5)
    pricers_smaller_sn_decompo = OrderedDict()
    smaller_sn_decompo_last_iter_sol = OrderedDict()
    smaller_sn_decompo_nb_desactivated = 0
    smaller_sn_decompo_overall_subpb = nb_part_smaller * length(vn_decompo.subgraphs)
    some_i = 1
    for vn_subgraph in vn_decompo.subgraphs
        pricers_smaller_sn_decompo[vn_subgraph] = set_up_pricer_sn_decompo(instance, vn_subgraph, smaller_sn_decompo_clusters)
        for pricer in pricers_smaller_sn_decompo[vn_subgraph]
            smaller_sn_decompo_last_iter_sol[pricer] = -99999. -some_i
            some_i += 1
        end
    end

    #----- sn partittioning bigger pricers
    nb_part_bigger = floor(Int, nv(s_network) / (nb_max_node_subgraphs*2.75))
    bigger_sn_decompo_clusters = get_sn_decompo(s_network, nb_part_bigger, nb_max_node_subgraphs*3.5)
    pricers_bigger_sn_decompo = OrderedDict()
    bigger_sn_decompo_last_iter_sol = OrderedDict()
    bigger_sn_decompo_nb_desactivated = 0
    bigger_sn_decompo_overall_subpb = nb_part_bigger * length(vn_decompo.subgraphs)
    some_i = 1
    for vn_subgraph in vn_decompo.subgraphs
        pricers_bigger_sn_decompo[vn_subgraph] = set_up_pricer_sn_decompo(instance, vn_subgraph, bigger_sn_decompo_clusters)
        for pricer in pricers_bigger_sn_decompo[vn_subgraph]
            bigger_sn_decompo_last_iter_sol[pricer] = -99999. -some_i
            some_i += 1
        end
    end


    #------ classical pricers
    pricers_full = Dict()
    for subgraph in vn_decompo.subgraphs
        pricers_full[subgraph] = set_up_pricer(instance, subgraph)
    end



    #------------ GENERATION DE COLONNES
    nb_columns = 0
    nb_iter = 0

    time_master = 0
    time_subproblems = 0
    time_overall = time()-time_beginning
    part_colge = 0.90
    time_colge = part_colge * time_max 
    time_end_solving = (1-part_colge) * time_max 

    keep_on = true
    algo_to_use = 1 # 1 is SN partitionning, 2 is reduced pricers, 3 is full pricers
    reason = "I don't know"
    CG_bound = 99999
    LG_bound = 0


    #----- Generating first columns => todo with heuristic in next round
    for vn_subgraph in vn_decompo.subgraphs
        for pricer_sub_sn in pricers_smaller_sn_decompo[vn_subgraph]
            update_pricer_sn_decompo(vn_decompo, pricer_sub_sn, dual_costs)
            column, obj_value = solve_pricers_sn_decompo(pricer_sub_sn)
            if obj_value < 9999
                add_column(master_problem, instance, vn_subgraph, column)
                nb_columns += 1
            end
        end
    end

    print("First set of $nb_columns columns generated, ")



    print("\n\n==================== Starting CG ====================\n")



    while keep_on
        nb_iter += 1
        print("Iter $nb_iter, ")

        optimize!(model)
        time_master +=  solve_time(model)

        status = termination_status(model)
        if status != MOI.OPTIMAL
            println("Infeasible or unfinished: $status")
            return
        end
        CG_bound = objective_value(model)
        print("CG value : $(round(CG_bound;   digits=3)), ")

        dual_costs = get_duals(instance, vn_decompo, master_problem)



        # ============ PRICER STUFF ============
        
        has_found_new_column = true
        time_beginning_pricer = time()

        # ------- SMALLER SN PARTITION ---------
        if algo_to_use == 1  
            print("pricer=smallerpaving, ")

            sorted_subpb = sort(collect(smaller_sn_decompo_last_iter_sol), by=x->x[2])
            overall_obj = 0


            for couple in sorted_subpb[1:5]
                pricer_sub_sn = couple[1]
                time_limit_subpb = time_colge - (time()-time_beginning)
                if time_limit_subpb < 0.01
                    break
                end

                update_pricer_sn_decompo(vn_decompo, pricer_sub_sn, dual_costs)
                column, obj_value = solve_pricers_sn_decompo(pricer_sub_sn, time_limit=time_limit_subpb)

                if column !== nothing
                    add_column(master_problem, instance, pricer_sub_sn.vn_subgraph, column)
                    nb_columns += 1
                end

                if obj_value>-5.
                    smaller_sn_decompo_nb_desactivated+=1
                end
                smaller_sn_decompo_last_iter_sol[pricer_sub_sn] = obj_value

                overall_obj += obj_value
            end
            print("mean sol: $(round(overall_obj/5;digits=3)), ")


            if smaller_sn_decompo_nb_desactivated >= smaller_sn_decompo_overall_subpb   
                print("too many subpb desactivated - changing pricers")
                algo_to_use=2
            end 


        elseif algo_to_use == 2  
            print("pricer=biggerpaving, ")

            sorted_subpb = sort(collect(bigger_sn_decompo_last_iter_sol), by=x->x[2])
            overall_obj = 0
            for couple in sorted_subpb[1:5]
                pricer_sub_sn = couple[1]

                time_limit_subpb = time_colge - (time()-time_beginning)
                if time_limit_subpb < 0.01
                    break
                end

                update_pricer_sn_decompo(vn_decompo, pricer_sub_sn, dual_costs)
                column, obj_value = solve_pricers_sn_decompo(pricer_sub_sn, time_limit=time_limit_subpb)

                if column !== nothing
                    add_column(master_problem, instance, pricer_sub_sn.vn_subgraph, column)
                    nb_columns += 1
                end

                if obj_value>-5.
                    bigger_sn_decompo_nb_desactivated+=1
                end
                bigger_sn_decompo_last_iter_sol[pricer_sub_sn] = obj_value

                overall_obj += obj_value
            end
            print("mean sol: $(round(overall_obj/5;digits=3)), ")


            if bigger_sn_decompo_nb_desactivated >= bigger_sn_decompo_overall_subpb   
                print("too many subpb desactivated - changing pricers")
                algo_to_use=3
            end 


                
        # --------- FULL PRICING -----------
        elseif algo_to_use == 3 # FULL PRICING
            print("exact, ")

            # do all subpb, get lg bound ?
            # also need to ensure that time is respected !
            has_found_new_column = false
            sum_pricers_values = 0
            for vn_subgraph in vn_decompo.subgraphs
                pricer = pricers_full[vn_subgraph]

                time_limit_subpb = time_colge - (time()-time_beginning)
                if time_limit_subpb < 0.01
                    break
                end

                column, obj_value = update_solve_pricer(instance, vn_decompo, pricer, dual_costs; time_limit = time_limit_subpb)

                if obj_value < -0.0001
                    has_found_new_column = true 
                    add_column(master_problem, instance, vn_subgraph, column)
                    nb_columns += 1
                    push!(vn_subgraph.columns, column)
                end

                sum_pricers_values += obj_value
            end

            # need to compute the LG bound here !
            LG_bound = CG_bound + sum_pricers_values
            print("LG bound: $(round(LG_bound; digits=3)), ")

        else
            println("pricer: not cooked yet! ")
            return
        end
        time_subproblems += time() - time_beginning_pricer

        print("with $nb_columns cols")



        time_overall = time()-time_beginning
        if time_overall < time_max * 0.90
            keep_on = true
            if algo_to_use == 2
                if !has_found_new_column
                    keep_on = false
                    reason="no improving columns"
                end
            end

        else
            keep_on = false
            reason="time limit"
        end

        print("in $(round(time_overall;digits=2))s")



        println("")
    end


        


    print("\n==================== CG finished ====================\nReason: $reason \n")
    println("Time in MP: $(round(time_master; digits=3)) , time in SP: $(round(time_subproblems; digits=3)), time overall: $(round(time_overall; digits=3))")
    println("$nb_iter iters, final value: $(round(CG_bound; digits=3))")
    println("====================================================\n")



    # ======= END HEURISTIC STUFF ======= #

    solution_heuristic = basic_heuristic(instance, vn_decompo, time_end_solving)


    return solution_heuristic, length(vn_decompo.subgraphs), CG_bound, LG_bound, nb_columns
    
end





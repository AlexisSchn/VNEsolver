# the best I can do for very large vn instances

using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX, Gurobi

#general
includet("../../utils/import_utils.jl")
includet("../../utils/visu.jl")

# utils colge
includet("utils/utils-subgraphdecompo.jl")
includet("utils/partition-vn.jl")

# pricers
includet("pricers/simplified-pricer.jl")
includet("pricers/pricer-full.jl")
includet("pricers/sn-decompo.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")



function solve_subgraph_decompo(instance, time_max = 20)

    println("Starting...")
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network


    # ======= SETTING UP THE DECOMPOSITION ======= #

    v_node_partitionning = partition_vn(instance)

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

    #----- sn partitionning pricers
    pricers_sn_decompo = Dict()
    for subgraph in vn_decompo.subgraphs
        pricers_sn_decompo[subgraph] = set_up_pricers_sn_decompo(instance, subgraph)
    end

    #------ reduced pricers
    pricer_reduced = Dict()
    for subgraph in vn_decompo.subgraphs
        pricer_reduced[subgraph] = set_up_pricer_reduced(instance, subgraph)
    end

    #------ classical pricers
    pricers_full = Dict()
    for subgraph in vn_decompo.subgraphs
        pricers_full[subgraph] = set_up_pricer(instance, subgraph)
    end






    #------------ GENERATION DE COLONNES
    nb_columns = 0
    time_master = 0
    time_subproblems = 0
    time_overall = time()-time_beginning
    keep_on = true
    algo_to_use = 1 # 1 is SN partitionning, 2 is reduced pricers, 3 is full pricers
    subgraph_to_do = 1 # Let's do them one by one since it's so simple...
    nb_iter = 0
    reason = "I don't know"
    CG_bound = 99999


    #----- Generating first columns
    for vn_subgraph in vn_decompo.subgraphs
        for pricer_sub_sn in pricers_sn_decompo[vn_subgraph]
            update_pricer_sn_decompo(vn_decompo, pricer_sub_sn, dual_costs)
            column, obj_value = solve_pricers_sn_decompo(pricer_sub_sn)
            if obj_value<999
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




        # ======= LAGRANGIAN BOUND STUFF =======
        
        
        





        # ======= PRICER STUFF =======

        vn_subgraph = vn_decompo.subgraphs[subgraph_to_do]
        print("pricing $(vn_subgraph.graph[][:name]) by ")
        
        
        has_found_new_column = false
        time_beginning_pricer = time()

        # ------- SN PARTITION ---------
        if algo_to_use == 1  

            best_obj_value = 9999
            print("sn paving, ")
            for pricer_sub_sn in pricers_sn_decompo[vn_subgraph]
                update_pricer_sn_decompo(vn_decompo, pricer_sub_sn, dual_costs)
                column, obj_value = solve_pricers_sn_decompo(pricer_sub_sn)
                if obj_value<-0.1
                    add_column(master_problem, instance, vn_subgraph, column)
                    has_found_new_column = true
                    nb_columns += 1
                    #println("obj value: $obj_value")
                    if obj_value < best_obj_value
                        best_obj_value = obj_value
                    end
                end
            end
            print("best sol: $(round(best_obj_value; digits=3)), ")


        # ---------- REDUCED PRICING ---------
        elseif algo_to_use == 2

            print("reduced, ")
            pricer = pricer_reduced[vn_subgraph]
            column, obj_value = update_solve_pricer_reduced(instance, vn_decompo, pricer, dual_costs)
            if obj_value < -0.0001
                has_found_new_column = true 
                add_column(master_problem, instance, vn_subgraph, column)
                nb_columns += 1
                push!(vn_subgraph.columns, column)
            end
            print("sol $(round(obj_value; digits=3)), ")




    
        # --------- FULL PRICING -----------
        elseif algo_to_use == 3 # FULL PRICING
            print("exact, ")
            pricer = pricers_full[vn_subgraph]
            column, obj_value = update_solve_pricer(instance, vn_decompo, pricer, dual_costs)
            if obj_value < -0.0001
                has_found_new_column = true 
                add_column(master_problem, instance, vn_subgraph, column)
                nb_columns += 1
                push!(vn_subgraph.columns, column)
            end
            print("sol $(round(obj_value; digits=3)), ")
        else
            println("pricer: not cooked yet! ")
            return
        end
        time_subproblems += time() - time_beginning_pricer

        print("with $nb_columns cols")


        # Next iteration cooking ?
        subgraph_to_do += 1
        if subgraph_to_do > length(vn_decompo.subgraphs)
            subgraph_to_do = 1
        end
        if nb_iter > 15
            algo_to_use = 2
        end
        if nb_iter > 20
            algo_to_use = 3
        end

        time_overall = time()-time_beginning
        if time_overall < time_max
            keep_on = true
        else
            keep_on = false
        end
        println("")
    end


        


    print("\n==================== CG finished ====================\nReason: $reason \n")
    println("Time in MP: $(round(time_master; digits=3)) , time in SP: $(round(time_subproblems; digits=3)), time overall: $(round(time_overall; digits=3))")
    println("$nb_iter iters, final value: $(round(CG_bound; digits=3))")
    println("====================================================\n")



    # ======= END HEURISTIC STUFF ======= #

    basic_heuristic(instance, vn_decompo)

    return
    
end



# todo for sauron
# => pave the network with cheapest heuristics
# => remove all the useless stuff, make the colge as clean as possible for ifip
# => better paving with cplex : get the right sizes that work better, get several solution at the same time


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
includet("pricers/basic-heuristic.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")



function solve_subgraph_decompo(instance; time_max = 100, v_node_partitionning = [], nb_part = -1)

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
            $(length(vn_decompo.v_edges_master)) cutting edges")

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

    # ====== STEP 1: simples heuristics paving / TO REPLACE WITH PAVING WITH ILP?
    println("First set of columns generated, or did I?")
    mappings = pave_heuristic(instance, vn_decompo)
    for v_subgraph in vn_decompo.subgraphs
        for truc in mappings[v_subgraph]
            add_column(master_problem, instance, v_subgraph, truc[1], truc[2])
        end
    end



    # ====== STEP 2: smaller SN pricers 
    nb_max_node_subgraphs = maximum([nv(subgraph.graph) for subgraph in vn_decompo.subgraphs])

    println("Ok so I have $(length(vn_decompo.subgraphs)), max has $nb_max_node_subgraphs, while sn has $(nv(s_network)).")
    nb_parts = floor(Int, nv(s_network) / (nb_max_node_subgraphs)/2)
    sn_decompo_clusters = get_sn_decompo(s_network, nb_parts, nb_max_node_subgraphs*2.5)
    println("And here is the nb cluster I ought to do: $nb_parts")
    pricers_sn_decompo = OrderedDict()
    sn_decompo_last_iter_sol = OrderedDict()
    nb_desactivated_pricers = 0
    nb_pricers = nb_parts * length(vn_decompo.subgraphs)
    for vn_subgraph in vn_decompo.subgraphs
        pricers_sn_decompo[vn_subgraph] = set_up_pricer_sn_decompo(instance, vn_subgraph, sn_decompo_clusters)
        for pricer in pricers_sn_decompo[vn_subgraph]
            sn_decompo_last_iter_sol[pricer] = -99999.
        end
    end

    keep_on = true
    reason = "I don't know"
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
        cg_value = objective_value(model)
        print("CG value : $(round(cg_value;   digits=3)), ")

        dual_costs = get_duals(instance, vn_decompo, master_problem)


        
        time_beginning_pricer = time()

        print("pricer=smaller-sn, ")

        sorted_subpb = sort(collect(sn_decompo_last_iter_sol), by=x->x[2])
        overall_obj = 0


        for couple in sorted_subpb[1:5]
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

            if reduced_cost>-5.
                nb_desactivated_pricers+=1
            end
            sn_decompo_last_iter_sol[pricer_sub_sn] = reduced_cost

            overall_obj += reduced_cost
        end
        print("mean sol: $(round(overall_obj/5;digits=3)), ")


 

        time_subproblems += time() - time_beginning_pricer

        print(", $nb_columns cols")



        time_overall = time()-time_beginning
        if time_overall < time_colge
            keep_on = true
            if nb_desactivated_pricers >= nb_pricers   
                keep_on = false
                reason="changing to full solving to get better columns"
            end 

        else
            keep_on = false
            reason="time limit"
        end

        print("in $(round(time_overall;digits=2))s")



        println("")
    end






    # ====== STEP 3: full pricers


    pricers_full = Dict()
    for subgraph in vn_decompo.subgraphs
        pricers_full[subgraph] = set_up_pricer(instance, subgraph)
    end


    keep_on = true
    reason = "I don't know"
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
        cg_value = objective_value(model)
        print("CG value : $(round(cg_value;   digits=3)), ")

        dual_costs = get_duals(instance, vn_decompo, master_problem)



        
        has_found_new_column = true
        time_beginning_pricer = time()

        print("exact, ")

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

        # need to compute the LG bound here !
        lower_bound = cg_value + sum_pricers_values
        print("LG bound: $(round(lower_bound; digits=3)), ")

        time_subproblems += time() - time_beginning_pricer

        print("with $nb_columns cols")



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

        print("in $(round(time_overall;digits=2))s")



        println("")
    end


        


    print("\n==================== CG finished ====================\nReason: $reason \n")
    println("Time in MP: $(round(time_master; digits=3)) , time in SP: $(round(time_subproblems; digits=3)), time overall: $(round(time_overall; digits=3))")
    println("$nb_iter iters, final value: $(round(cg_value; digits=3))")
    println("====================================================\n")



    # ======= END HEURISTIC STUFF ======= #

    solution_heuristic = basic_heuristic(instance, vn_decompo, master_problem, time_end_solving)

    return solution_heuristic
    
end



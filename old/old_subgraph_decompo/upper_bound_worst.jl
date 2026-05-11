using Printf
using Graphs, MetaGraphsNext
using JuMP, CPLEX


#general
includet("../../utils/import_utils.jl")

# utils colge
includet("utils/master-problem.jl")
includet("utils/graph-decomposition.jl")
includet("utils/partition-graph.jl")

# pricers
includet("pricers/milp-pricer.jl")
includet("pricers/local-search-pricer.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")



function solve_price_branch_bad(instance; pricer="milp", nb_virtual_subgraph=0, nb_columns_max=300, time_end_milp=30)

    # === SOME USEFUL THINGS === #
    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    println("Starting...")
    time_beginning = time()
    nb_columns = 0


    # === SOME PARAMETERS === #
    alpha_colge=0.5
    nb_columns_init = min(nb_columns_max, 50)




    # ======= SETTING UP THE DECOMPOSITION ======= #

    # ------ Virtual decomposition ------ #
    if nb_virtual_subgraph == 0
        nb_virtual_subgraph = floor(Int, nv(v_network.graph)/10)
    end
    v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.2)   
    vn_decompo = set_up_decompo(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs
    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")
    



    # Starting the algorithm
    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network_dir)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir, s_edge, :cost)
    end
    shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)

    empty_dual_costs = get_empty_duals(instance, vn_decompo)


    # ======= INITIALIZATION  ======= #
    println("STEP 1: INITIALIZATION...")
    time_0=time()

    result = find_columns_dumb(  instance, vn_subgraphs, vn_decompo, empty_dual_costs;  solver=pricer)
    sub_mappings = result[:mappings]
    for v_subgraph in vn_decompo.subgraphs
        for mapping in sub_mappings[v_subgraph]
            add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
            nb_columns += 1
        end
    end

    println("\nInitialization dong: \n$nb_columns columns gotten, with $pricer, in $(time() - time_0)s\n")


    # ======= FIRST ITERATION OF MASTER PROBLEM ======= #
    optimize!(model)
    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("Master problem is infeasible or unfinished: $status")
        return
    end
    rmp_value = objective_value(model)
    time_overall = time() - time_beginning
    nb_iter = 0
    @printf("Iter %2d  RMP value: %10.3f  %5d column    time: %5.2fs  \n", nb_iter, rmp_value, nb_columns, time_overall)


    # ======= COLUMN GENERATION ======= #
    while nb_columns < nb_columns_max

        # --- pricers
        current_dual_costs = get_duals(instance, vn_decompo, master_problem)
        dual_costs_to_use = average_dual_costs(instance, vn_decompo, empty_dual_costs, current_dual_costs, alpha=alpha_colge)
        
        result =  find_columns_dumb(  instance, vn_subgraphs, vn_decompo, dual_costs_to_use;  solver=pricer)
        sub_mappings = result[:mappings]
        for v_subgraph in vn_decompo.subgraphs
            for mapping in sub_mappings[v_subgraph]
                add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
                nb_columns += 1
            end
        end


        # ---- master problem part
        optimize!(model)
        rmp_value = objective_value(model)


        # --- printing and utils
        time_overall = time()-time_beginning
        nb_iter += 1
        #average_obj = overall_reduced_costs / (nb_feasible+10e-6)
        #=
        @printf("Iter %2d    RMP value:%10.3f  %5d columns  time: %5.2fs  mean red. cost: %10.3f   feasible %10.3f / %10.3f \n",
            nb_iter, rmp_value, nb_columns, time_overall, average_obj, nb_feasible, nb_virtual_subgraph
        )
        =#
        @printf("Iter %2d    RMP value:%10.3f  %5d columns  time: %5.2fs \n",
            nb_iter, rmp_value, nb_columns, time_overall
        )
    end

    # ======= END HEURISTICS ======= #

    # ---- Restricted master problem heuristic
    value_cg_heuristic, cg_heuristic_solution = basic_heuristic(instance, vn_decompo, master_problem, time_end_milp)

    println("\nSolution obtained: $value_cg_heuristic\n")

    #local_search_changin(instance, cg_heuristic_solution, 300)

    return (mapping_cost=value_cg_heuristic,
        mapping = nothing  )
end







function find_columns_dumb(  instance, vn_subgraphs, vn_decompo, dual_costs;  solver="milp")


    # Useful things
    mappings_per_subgraph = Dict()

    for v_subgraph in vn_subgraphs
        mappings_per_subgraph[v_subgraph] = []
        # GETTING THE SUBMAPPING
        if solver == "local-search"
            sub_instance = Instance(v_subgraph.graph, instance.s_network)
            result =    solve_local_search_pricer(sub_instance, dual_costs; nb_iterations=250)
            sub_mapping = result[:sub_mapping]
            cost = result[:real_cost]
        elseif solver == "milp"
            pricer = set_up_pricer(instance, v_subgraph)
            sub_mapping, cost, reduced_costs = update_solve_pricer(instance, vn_decompo, pricer, dual_costs)
            result = update_solve_pricer(instance, vn_decompo, pricer, dual_costs, time_limit = 10)
        else
            println("Pricer unknown???")
            return
        end

        if isnothing(sub_mapping) # invalid submapping!
            #print("A pricer failed. ")
            continue
        end
    
        push!(mappings_per_subgraph[v_subgraph], sub_mapping)
        

    end


    return (    mappings=mappings_per_subgraph, 
                average_dual_costs=average_dual_costs
    )    


end



using Printf
using Graphs, MetaGraphsNext
using JuMP, CPLEX


#general
includet("../../utils/import_utils.jl")

# utils colge
includet("utils/master-problem.jl")
includet("utils/graph-decomposition.jl")
includet("utils/partition-graph.jl")

# init
includet("init/init-paving-routing.jl")

# pricers
includet("pricers/milp/pricer-partition-sn-routing-penalty.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")
includet("end-heuristic/local-search-exact.jl")



function solve_drake(instance; nb_virtual_subgraph=0, nb_columns_max=300)



    # Budget : 60 seconds
    time_submappings = 30
    time_cg_heuristic = 30

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    println("Starting...")
    time_beginning = time()
    nb_columns = 0


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
    
    # ------ Substrate decomposition ------ #
    size_max_v_subgraph = maximum(nv(v_subgraph.graph) for v_subgraph in vn_decompo.subgraphs)
    nb_substrate_nodes_capacited = 0
    for s_node in vertices(s_network)
        if get_attribute_node(s_network, s_node, :cap) >= 1
            nb_substrate_nodes_capacited += 1
        end
    end
    ratio_capacited = nb_substrate_nodes_capacited / nv(s_network)
    nb_substrate_subgraphs = maximum([floor(Int, nv(s_network) / (size_max_v_subgraph*1.5/ratio_capacited)), nb_virtual_subgraph])
    clusters = partition_graph(s_network.graph, nb_substrate_subgraphs; max_umbalance = 1.3)
    sn_subgraphs = []
    for (i_subgraph, cluster) in enumerate(clusters)
        induced_subg = my_induced_subgraph(s_network, cluster, "sub_sn_$i_subgraph")
        push!(sn_subgraphs, Subgraph(induced_subg, cluster))
    end
    println("Substrate network decomposition done:")
    print_stuff_subgraphs(s_network, sn_subgraphs)



    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model



    # ======= PAVING THE NETWORK  ======= #
    println("STEP 1: INITIALIZATION...")
    nb_mappings_first = 100

    nb_cols_notnew = 0

    time_0=time()
    sub_mappings = find_submappings_routing(instance, vn_decompo, sn_subgraphs, solver="milp", nb_columns=nb_mappings_first)

    for v_subgraph in vn_decompo.subgraphs
        for mapping in sub_mappings[v_subgraph]
            if !check_if_column_new(master_problem, mapping, v_subgraph)
                nb_cols_notnew += 1
            end
            add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
            nb_columns += 1
        end
    end

    println("$nb_columns mappings gotten! In just $(time() - time_0)")


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

    @printf("Iter %2d  RMP value: %10.3f  %5d column    time: %5.2fs  \n",
        nb_iter, rmp_value, nb_columns, time_overall
    )

    


    # ======= PAVING THE NETWORK WITH COLUMN GENERATION ======= #

    # Base shortest paths
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network_dir)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir, s_edge, :cost)
    end
    base_shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)

    alpha_colge=0.
    used_dual_costs = get_empty_duals(instance, vn_decompo)
    pricers_sn_partition = set_up_pricers_sn_partitionning(instance, vn_subgraphs, sn_subgraphs)

    while nb_columns < nb_columns_max

        # ---- pricers part

        # Associate subvn to a random subsn
        used_sub_s_network = []
        assignment_virtual_substrate_subgraphs = Dict()
        for v_subgraph in vn_subgraphs
            found = false
            while !found
                i_subgraph = rand(1:nb_substrate_subgraphs)
                if i_subgraph ∉ used_sub_s_network
                    push!(used_sub_s_network, i_subgraph)
                    found = true
                    assignment_virtual_substrate_subgraphs[v_subgraph] = sn_subgraphs[i_subgraph]
                end
            end
        end


        # Initialize temporary placement
        temporary_placement = zeros(Int, nv(v_network))
        for v_subgraph in vn_subgraphs
            s_subgraph = assignment_virtual_substrate_subgraphs[v_subgraph]
            cluster = s_subgraph.nodes_of_main_graph
            capacities_score = [ get_attribute_node(s_network, s_node, :cap) * 
                    sum(get_attribute_edge(s_network, get_edge(s_network, s_node, s_neighbor), :cap) for s_neighbor in neighbors(s_network, s_node)) 
                    for s_node in cluster]
            best_node = cluster[findmin(capacities_score)[2]]
            for v_node in v_subgraph.nodes_of_main_graph
                temporary_placement[v_node] = best_node
            end
        end
    
        
        # Solve
        old_dual_costs = get_empty_duals(instance, vn_decompo)
        #old_dual_costs = used_dual_costs
        current_dual_costs = get_duals(instance, vn_decompo, master_problem)
        used_dual_costs = average_dual_costs(instance, vn_decompo, old_dual_costs, current_dual_costs, alpha=alpha_colge)
        
        nb_feasible = 0
        overall_reduced_costs = 0
        for v_subgraph in vn_subgraphs
            s_subgraph = assignment_virtual_substrate_subgraphs[v_subgraph]


            # Additional routing costs
            additional_costs = []
            for v_node in vertices(v_subgraph.graph)
                current_addition_costs = [0 for s_node in vertices(s_subgraph.graph)]
                original_v_node = v_subgraph.nodes_of_main_graph[v_node]
                
                for v_edge in vn_decompo.v_edges_master
                    if src(v_edge) == original_v_node
                        placement_of_dst_node = temporary_placement[dst(v_edge)]
                        for s_node in vertices(s_subgraph.graph)
                            original_s_node = s_subgraph.nodes_of_main_graph[s_node]
                            current_addition_costs[s_node] += base_shortest_paths.dists[original_s_node, placement_of_dst_node] 
                        end
                    end
                    if dst(v_edge) == original_v_node
                        placement_of_dst_node = temporary_placement[src(v_edge)]
                        for s_node in vertices(s_subgraph.graph)
                            original_s_node = s_subgraph.nodes_of_main_graph[s_node]
                            current_addition_costs[s_node] += base_shortest_paths.dists[original_s_node, placement_of_dst_node] 
                        end
                    end
                end
                
    
                push!(additional_costs, current_addition_costs)
            end

                        
            pricer = pricers_sn_partition[v_subgraph][s_subgraph]
            update_pricer_sn_decompo_penalty(vn_decompo, pricer, used_dual_costs, additional_costs)
            result_pricer = solve_pricers_sn_decompo_penalty(pricer, additional_costs)
            if result_pricer["true_cost"] < 10e6
                
                nb_feasible+=1
                overall_reduced_costs += result_pricer["reduced_cost"]
                sub_mapping = result_pricer["mapping"]
                if !check_if_column_new(master_problem, sub_mapping, v_subgraph)
                    nb_cols_notnew += 1
                end
                add_column(master_problem, instance, v_subgraph, sub_mapping, result_pricer["true_cost"])
                nb_columns += 1

                if (nb_columns % 100) == 0
                    println("\nSo far: $nb_cols_notnew / $nb_columns\n")
                end

                # Temporary placement stuff
                for v_node in vertices(v_subgraph.graph)
                    real_s_node = sub_mapping.node_placement[v_node]
                    real_v_node = v_subgraph.nodes_of_main_graph[v_node]
                    temporary_placement[real_v_node] = real_s_node
                end

            end
        end





        # ---- master problem part
        optimize!(model)
        rmp_value = objective_value(model)


        # --- printing and utils
        time_overall = time()-time_beginning
        nb_iter += 1
        average_obj = overall_reduced_costs / (nb_feasible+10e-6)

        @printf("Iter %2d    RMP value:%10.3f  %5d columns  time: %5.2fs  mean red. cost: %10.3f   feasible %10.3f / %10.3f \n",
            nb_iter, rmp_value, nb_columns, time_overall, average_obj, nb_feasible, nb_virtual_subgraph
        )

    end

    # ======= END HEURISTICS ======= #

    # ---- Restricted master problem heuristic
    time_cg_heuristic = 60
    value_cg_heuristic, cg_heuristic_solution = basic_heuristic(instance, vn_decompo, master_problem, time_cg_heuristic)

    println("\nSolution obtained: $value_cg_heuristic\n")

    #local_search_changin(instance, cg_heuristic_solution, 300)

    return Dict("mapping_cost"=>value_cg_heuristic)

end









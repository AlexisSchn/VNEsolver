using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf


#general
includet("../../utils/import_utils.jl")

# utils colge
includet("utils/utils-subgraphdecompo.jl")
includet("utils/partition-graph.jl")
includet("utils/checkers.jl")


# init
includet("init/init-paving-routing.jl")
includet("init/init-paving-simple.jl")

# pricers
includet("pricers/milp/pricer-subsn.jl")
includet("pricers/milp/pricer-exact.jl")
includet("pricers/mepso/mepso-pricer-subsn.jl")
includet("pricers/local-search/ls-pricer-subsn.jl")


function crazy_bound(instance)
    
    
    println("Starting...")

    time_master = 0
    nb_iter=0
    lower_bound = 0
    nb_columns=0

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    time_beginning = time()


    # ======= SETTING UP THE DECOMPOSITION ======= #
    nb_virtual_subgraph = floor(Int, nv(v_network.graph)/10)
    v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.2)

    vn_decompo = set_up_decompo(instance, v_node_partitionning)
    
    println("Decomposition set: ")
        println("For $v_network, there is $(length(vn_decompo.subgraphs)) subgraphs:")

    for subgraph in vn_decompo.subgraphs
        println("       $(subgraph.graph[][:name]) with $(nv(subgraph.graph)) nodes, $(ne(subgraph.graph)) edges")
    end
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")

    
    vn_subgraphs = vn_decompo.subgraphs

    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    print("Master problem set... ")

    # ====== STEP 1 : INIT ======= #
    nb_columns_to_add_init = 200

    println("Paving time...")
    time_0 = time()

    # Get substrate subgraphs - 1
    size_max_v_subgraph = maximum(nv(v_subgraph.graph) for v_subgraph in vn_decompo.subgraphs)
    nb_substrate_subgraphs = floor(Int, nv(s_network) / (size_max_v_subgraph*1.5))
    nb_nodes_subgraph = 20
    #clusters = get_sn_decompo(s_network, nb_substrate_subgraphs, nb_nodes_subgraph)

    clusters = partition_graph(s_network.graph, nb_substrate_subgraphs; max_umbalance = 1.3)
    sn_subgraphs = []
    for (i_subgraph, cluster) in enumerate(clusters)
        print("Cluster $i_subgraph has $(length(cluster)) nodes ")
        induced_subg = my_induced_subgraph(s_network, cluster, "sub_sn_$i_subgraph")
        push!(sn_subgraphs,Subgraph(induced_subg, cluster))
    end

    #sub_mappings = find_submappings_routing(instance, vn_decompo, sn_subgraphs, solver="mepso", nb_columns=nb_columns_to_add_init)
    sub_mappings = find_submappings_simple(instance, vn_decompo, sn_subgraphs, solver="local-search", nb_columns=nb_columns_to_add_init)

    println("Mappings gotten! In just $(time() - time_0)")
    for v_subgraph in vn_decompo.subgraphs
        for mapping in sub_mappings[v_subgraph]
            add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
            nb_columns+=1
        end
    end
    print("Submappings added - 1...")



    # ======= ITERATION OF MASTER PROBLEM ======= #
    optimize!(model)
    time_master +=  solve_time(model)

    cg_value = objective_value(model)
    time_overall = time() - time_beginning
    @printf("Iter %2d  RMP value: %10.3f  %5d column    time: %5.2fs  \n",
        nb_iter, cg_value, nb_columns, time_overall
    )


    # ======= FIRST ITERATION OF MASTER PROBLEM ======= #
    optimize!(model)
    time_master +=  solve_time(model)

    cg_value = objective_value(model)
    time_overall = time() - time_beginning
    @printf("Iter %2d  RMP value: %10.3f  %5d column    time: %5.2fs  \n",
        nb_iter, cg_value, nb_columns, time_overall
    )


    # ====== STEP 2: smaller pricers - paving the network
    println("------- Part 2: Reduced pricers")




    println("YAY I MEAN LET'S LOOK AT THE DUAL COSTS EHEHEHEHHEHEHEHEHEHEH")
    dual_costs = get_duals(instance, vn_decompo, master_problem)
    
    capacity_s_node_history = []
    for s_node in vertices(s_network)
        push!(capacity_s_node_history, [dual_costs.capacity_s_node[s_node]])
    end

    flow_conservation_history= Dict()
    for v_edge in vn_decompo.v_edges_master
        
        flow_conservation_history[v_edge] = []
        for s_node in vertices(s_network)
            push!(flow_conservation_history[v_edge], [dual_costs.flow_conservation[v_edge][s_node]])
        end
    end

    departure_history= Dict()
    for v_edge in vn_decompo.v_edges_master
        departure_history[v_edge] = []
        for s_node in vertices(s_network)
            push!(departure_history[v_edge], [dual_costs.departure[v_edge][s_node]])
        end
    end



    nb_columns_to_put = 300
    nb_columns_smallpricer = 0

    nb_substrate_subgraph = floor(Int, nv(s_network) / 15)  
    nb_nodes_subgraph = 27
    sn_decompo_clusters = get_sn_decompo(s_network, nb_substrate_subgraph, nb_nodes_subgraph)
    println("We have $nb_substrate_subgraph sub-substrate, with at least $nb_nodes_subgraph capacited nodes")

    pricers_sn_decompo = OrderedDict()
    sub_pricers_last_values = OrderedDict()
    for vn_subgraph in vn_decompo.subgraphs
        pricers_sn_decompo[vn_subgraph] = set_up_pricer_sn_decompo(instance, vn_subgraph, sn_decompo_clusters)
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

        # IM WATCHINGGGG

        for s_node in vertices(s_network)
            push!(capacity_s_node_history[s_node], dual_costs.capacity_s_node[s_node])
        end
        
        for v_edge in vn_decompo.v_edges_master
            for s_node in vertices(s_network)
                push!(flow_conservation_history[v_edge][s_node], dual_costs.flow_conservation[v_edge][s_node])
                push!(departure_history[v_edge][s_node], dual_costs.departure[v_edge][s_node])
            end
        end
    
        
        for key in keys(sub_pricers_last_values)
            sub_pricers_last_values[key] = sub_pricers_last_values[key]*1.1
        end
        
        sorted_subpb = sort(collect(sub_pricers_last_values), by=x->x[2])
        average_obj = 0
        nb_pricer_to_do = min(5, length(keys(sub_pricers_last_values)))
        for couple in sorted_subpb[1:nb_pricer_to_do]


            pricer_sub_sn = couple[1]


            update_pricer_sn_decompo(vn_decompo, pricer_sub_sn, dual_costs)
            column, true_cost, reduced_cost = solve_pricers_sn_decompo(pricer_sub_sn, time_limit=50)

            if column !== nothing && reduced_cost < -0.001
                add_column(master_problem, instance, pricer_sub_sn.vn_subgraph, column, true_cost)
                nb_columns += 1
                nb_columns_smallpricer += 1
            end

            sub_pricers_last_values[pricer_sub_sn] = reduced_cost

            average_obj += (reduced_cost/nb_pricer_to_do)
        end


        # ---- master problem part

        optimize!(model)
        time_master +=  solve_time(model)       
        println("Time last iter master : $(solve_time(model))")
        cg_value = objective_value(model)

        time_overall = time()-time_beginning

        @printf("Iter %2d  CG bound: %10.3f  lower bound: %10.3f  %5d column  time: %5.2fs  average reduced cost: %10.3f \n",
                    nb_iter, cg_value, lower_bound, nb_columns, time_overall, average_obj)
    


        
        # ----- useful things

        if nb_columns_smallpricer > nb_columns_to_put
            keep_on = false
        end

    end
    println("\n Step 2 finished, reason: $reason.")
    


    println("AND NOW MADAME MONSIEUR:")


    println("Capacity node history:")
    for s_node in vertices(s_network)
        println("Node $s_node: $(capacity_s_node_history[s_node])")
    end



    println("Flow conservation history:")
    for v_edge in vn_decompo.v_edges_master
        println("\n FOR $v_edge \n\n")
        for s_node in vertices(s_network)
            #println("for node $s_node: $(flow_conservation_history[v_edge][s_node])")
            println("For node $s_node: minimum $(minimum(flow_conservation_history[v_edge][s_node])), maximum $(maximum(flow_conservation_history[v_edge][s_node]))")
        end
    end


    println("\n\n\n\n\nDeparture history:")
    for v_edge in vn_decompo.v_edges_master
        println("\n FOR $v_edge \n\n")
        for s_node in vertices(s_network)
            #println("for node $s_node: $(departure_history[v_edge][s_node])")
            println("For node $s_node: minimum $(minimum(departure_history[v_edge][s_node])), maximum $(maximum(departure_history[v_edge][s_node]))")
        end
    end


    print("\n==================== CG finished ====================\nReason: $reason \n")
    println("Time in MP: $(round(time_master; digits=3)) , time in SP: idk, time overall: $(round(time_overall; digits=3))")
    println("$nb_iter iters, final value: $(round(cg_value; digits=3))")
    println("====================================================\n")



    
    
    # ======= GETTING A SOLUTION ======= #
    # LETS LOOK AT THE THINGS WE HAVE EHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEHEH

    println("Looking at the solution:")

    for v_subgraph in vn_decompo.subgraphs
        println("For $v_subgraph, we have:")
        for column in master_problem.columns[v_subgraph]
            if value(column.variable) > 0.0001
                println("Column $column with $(value(column.variable))")
            end
        end
        println("\n\n\n")
    end
    

    result = Dict()
    result["solving_time"] = time() - time_beginning
    result["mapping_cost"] = value_cg_heuristic

    return result
end







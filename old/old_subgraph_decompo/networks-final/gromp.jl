using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf


#general
includet("../../utils/import_utils.jl")

# utils colge
includet("utils/master-problem.jl")
includet("utils/graph-decomposition.jl")
includet("utils/partition-graph.jl")

# heuristics
includet("init/init-paving-routing.jl")
includet("init/init-paving-simple.jl")

# Find a final solution
includet("end-heuristic/basic-ilp.jl")
includet("end-heuristic/local-search-exact.jl")



function solve_gromp(instance; nb_virtual_subgraph=0, nb_submappings=150, routing_penalty=true)
    
    v_network = instance.v_network
    s_network = instance.s_network

    println("Starting...")
    time_beginning = time()



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







    # ====== PAVING THE NETWORK WITH HEURISTIC ======= #

    println("Paving time...")
    time_0 = time()

    if routing_penalty
        sub_mappings = find_submappings_routing(instance, vn_decompo, sn_subgraphs, solver="local-search", nb_columns=nb_submappings)
    else
        sub_mappings = find_submappings_simple(instance, vn_decompo, sn_subgraphs, solver="local-search", nb_columns=nb_submappings)
    end

    println("Mappings gotten! In just $(time() - time_0)")
    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    print("Master problem set... ")
    for v_subgraph in vn_decompo.subgraphs
        for mapping in sub_mappings[v_subgraph]
            check_if_column_new(master_problem, mapping, v_subgraph)
            add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
        end
        println("Finished with that subgraph...")
    end
    print("Submappings added...")

    

    # ======= GETTING A SOLUTION ======= #

    # For info, i print the bound
    optimize!(model)
    rmp_value = objective_value(model)
    println("RMP value: $(rmp_value)")
    time_cg_heuristic = 60
    value_cg_heuristic, cg_heuristic_solution = basic_heuristic(instance, vn_decompo, master_problem, time_cg_heuristic)

    return Dict(
        "solving_time" => (time() - time_beginning),
        "mapping_cost" => value_cg_heuristic,
        "rmp_value" => rmp_value
    )
    
end
















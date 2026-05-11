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

# pricers
includet("init/init_uepso.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")
includet("end-heuristic/local-search-exact.jl")



function solve_chill(instance)
    
    
    # Budget : 60 seconds
    time_init = 30
    time_cg_heuristic = 30


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




    # ====== PAVING THE NETWORK WITH HEURISTIC ======= #

    println("Paving time...")
    time_0 = time()
    # max 300 columns, or 50 per subgraphs.
    nb_column_per_subgraph = floor(Int,400/nb_virtual_subgraph)
    mappings = init_uepso(instance, vn_decompo, nb_column_per_subgraph)
    println("Mappings gotten! In just $(time() - time_0)")
    for v_subgraph in vn_decompo.subgraphs
        for mapping in mappings[v_subgraph]
            add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
        end
    end


        

    
    # ======= GETTING A SOLUTION ======= #
    value_cg_heuristic, cg_heuristic_solution = basic_heuristic(instance, vn_decompo, master_problem, time_cg_heuristic)
    #local_search(instance, vn_decompo, heur_sol)



    result = Dict()
    result["algo"] = "chiller-subgraph-decompo"
    result["time_solving"] = time() - time_beginning
    result["value_cg_heuristic"] = value_cg_heuristic

    return result
end



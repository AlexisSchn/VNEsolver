# A clean version...

using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf

#general
includet("../utils/import_utils.jl")

# utils colge
includet("utils/master_problem.jl")
includet("utils/graph_decomposition_overlapping.jl")
includet("utils/column_generation.jl")
includet("../utils/partition-graph.jl")

# init
includet("init/greedy_easy.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")



function strict_partition(instance; v_node_partitionning = [], nb_virtual_subgraph=0)

    println("Starting...")
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    # ======= SETTING UP THE DECOMPOSITION ======= #

    # AUTOMATIC PARTITION
    if v_node_partitionning == []
        if nb_virtual_subgraph == 0
            nb_virtual_subgraph = floor(Int, nv(v_network.graph)/10)
        end
        v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.2)
    end

    println("Node partitionning: $v_node_partitionning")

    vn_decompo = set_up_decompo_overlapping(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    
    # === COLUMN GENERATION === #

    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")



    # generating first columns. Completly random...
    nb_columns = 0
    while nb_columns < 300
        given_placement = Dict()
        for vn_subgraph in vn_subgraphs
            sub_instance = Instance(vn_subgraph.graph, s_network, s_network_dir)
            result = solve_greedy(sub_instance, given_placement)
            cost = result[:mapping_cost]
            if cost < 10e6
                sub_mapping = result[:mapping]
                add_column(master_problem, instance, vn_decompo, vn_subgraph, sub_mapping, cost)
                nb_columns +=1
            end
        end
    end
    println("$nb_columns generated before the CG!")

    # column generation!
    column_generation(instance, vn_decompo, master_problem)


    # ======= END HEURISTIC STUFF ======= #

    basic_heuristic(instance, vn_decompo, master_problem, 100)

    return 
end



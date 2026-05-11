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



function automatic_overlapping(instance; v_node_partitionning = [], nb_virtual_subgraph=0, nb_overlapping_nodes=3)

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
        v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.3)
    end

    println("Node partitionning: $v_node_partitionning")

    subgraphs = []
    for part in v_node_partitionning
        sg, mapping = induced_subgraph(v_network, part)
        edges_of_subg = collect(edges(sg))
        real_edges = []
        for edge in edges_of_subg
            push!(real_edges, get_edge(v_network, mapping[src(edge)], mapping[dst(edge)]))
        end
        println("Edges : $edges_of_subg, real edges $real_edges, mapping $mapping")
        push!(subgraphs, Dict("nodes"=>part, "edges"=>real_edges))
    end


    # Copy of the v network, kinda useful

    copy_v_network = copy(v_network.graph)

    for subg in subgraphs
        for v_edge in subg["edges"]
            rem_edge!(copy_v_network, v_edge)
        end
    end


    # Here, add the overlapping nodes!
    nb_nodes_added = 0
    

    while nb_nodes_added < nb_overlapping_nodes
        best_node = 0
        best_subgraph = 0
        best_score = 0
        
        for v_node in vertices(v_network)
            for (i_cluster, cluster) in enumerate(v_node_partitionning)
                if v_node ∉ cluster
                    nb_neighbor = length(neighbors(copy_v_network, v_node) ∩ v_node_partitionning[i_cluster] )
                    if nb_neighbor > best_score
                        best_score = nb_neighbor
                        best_node = v_node
                        best_subgraph = i_cluster
                    end
                end
            end
        end

        if best_score > 1
            println("I'm adding $best_node to cluster $best_subgraph, it had $best_score neighbors there!")
            println("Also, time to remove some edges...")
            for v_node in v_node_partitionning[best_subgraph]
                if v_node ∈ neighbors(copy_v_network, best_node)
                    edge = get_edge(v_network, v_node, best_node)
                    print(" another one found! $edge")
                    push!(subgraphs[best_subgraph]["edges"], edge)
                    rem_edge!(copy_v_network, edge)
                end
            end
            #push!(v_node_partitionning[best_subgraph], best_node)
            push!(subgraphs[best_subgraph]["nodes"], best_node)
            nb_nodes_added += 1
        else
            println("Didnt find any node... weird...")
            break
        end
    end


    println(" And in the end... the subgraphs : $subgraphs")

    vn_decompo = set_up_decompo_overlapping_more_info(instance, subgraphs)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    
    # === COLUMN GENERATION === #

    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")




    # column generation!
    return column_generation(instance, vn_decompo, master_problem)


    # ======= END HEURISTIC STUFF ======= #

    #basic_heuristic(instance, vn_decompo, master_problem, 900)

    return 
end



function automatic_overlapping_old(instance; v_node_partitionning = [], nb_virtual_subgraph=0, nb_overlapping_nodes=3)

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
        v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.3)
    end

    println("Node partitionning: $v_node_partitionning")

    # Here, add the overlapping nodes!
    nb_nodes_added = 0

    original_clusters = deepcopy(v_node_partitionning) # It will be simpler this way... 

    while nb_nodes_added < nb_overlapping_nodes
        best_node = 0
        best_subgraph = 0
        best_score = 0
        
        for v_node in vertices(v_network)
            for (i_cluster, cluster) in enumerate(v_node_partitionning)
                if v_node ∉ cluster
                    nb_neighbor = length(neighbors(v_network, v_node) ∩ original_clusters[i_cluster] )
                    if nb_neighbor > best_score
                        best_score = nb_neighbor
                        best_node = v_node
                        best_subgraph = i_cluster
                    end
                end
            end
        end

        if best_score > 1
            println("I'm adding $best_node to cluster $best_subgraph, it had $best_score neighbors there!")
            push!(v_node_partitionning[best_subgraph], best_node)
            nb_nodes_added += 1
        else
            println("Didnt find any node... weird...")
            break
        end
    end



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



    # generating first columns. Adapted for overlapping cg...
    #=
    node_capacities = [get_attribute_node(s_network, s_node, :cap) for s_node in vertices(s_network)]
    capacited_nodes = [s_node for s_node in vertices(s_network) if node_capacities[s_node] ≥ 1]
    centrality_nodes = closeness_centrality(s_network)
    nb_columns = 0
    while nb_columns < 300
        # Generate a random mapping for each subvn
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

        # Also, generate a random location for each overlapping node (on a pretty central substrate node)
        # This should help convergence?
        # The other method would be to add to the previous one somethin similar..
        for overlapping_node in keys(vn_decompo.overlapping_nodes)
            s_nodes_scores = [ (centrality_nodes[s_node] + 0.5 * rand() ) for s_node in capacited_nodes ]
            s_node = capacited_nodes[argmin(s_nodes_scores)]
            for vn_subgraph in keys(vn_decompo.v_nodes_assignment[overlapping_node])
                given_placement = Dict()
                v_node_in_subgraph = vn_decompo.v_nodes_assignment[overlapping_node][vn_subgraph]
                given_placement[v_node_in_subgraph] = s_node
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
    end
    println("$nb_columns generated before the CG!")
    =#

    # column generation!
    return column_generation(instance, vn_decompo, master_problem)


    # ======= END HEURISTIC STUFF ======= #

    #basic_heuristic(instance, vn_decompo, master_problem, 900)

    return 
end

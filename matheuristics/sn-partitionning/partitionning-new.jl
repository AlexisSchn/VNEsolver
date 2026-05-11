using Base.Libc


using Revise, JuMP, CPLEX, Gurobi, Metis
using Graphs, MetaGraphsNext

includet("../../../utils/import_utils.jl")
includet("../../../utils/visu.jl")
includet("kahip_wrapper.jl")
includet("../../compact_undir.jl")





function solve_with_partitionning(instance, time_limit = 3)

    time_beginning = time()

    s_network = instance.s_network
    v_network = instance.v_network

    # 1 : Partitionner
    nb_clusters = floor(Int, nv(s_network)/nv(v_network)) + 1
    inbalance = 0.
    println("$nb_clusters clusters to do, with inbalance $inbalance...")
    partition = partition_kahip(s_network.graph, nb_clusters, inbalance)
    clusters = [Vector{Int64}() for i in 1:nb_clusters]
    for s_node in vertices(s_network)
        push!(clusters[partition[s_node]], s_node)
    end

    # 2 : Corriger

    # 2a) correction by removing unconnected sets
    for cluster in clusters
        simple_subgraph, vmap = induced_subgraph(s_network.graph, cluster)
        if !is_connected(simple_subgraph)
            #println("It's not connected ! :(")
            #println("At the beginning the cluster is : $cluster")
            components = connected_components(simple_subgraph)
            component_sorted = sort(components, by=x->length(x), rev=true)
            new_cluster = [vmap[i] for i in component_sorted[1]]
            for subcluster in component_sorted[2:length(component_sorted)]
                #Let's add all those nodes to a (most) connected subgraph
                nodes_original = [vmap[node] for node in subcluster]
                subgraph_neighbors = zeros(Int, nb_clusters)
                for node in nodes_original
                    for neighbor in neighbors(s_network, node)
                        if neighbor ∉ cluster
                            subgraph_neighbors[partition[neighbor]] += 1
                        end
                    end
                end
                most_connected_subgraph = sortperm(subgraph_neighbors, rev=true)
                append!(clusters[most_connected_subgraph[1]], nodes_original)
                #println("Well let's add $nodes_original to cluster $(clusters[most_connected_subgraph[1]])")
                filter!(e->e∉nodes_original, cluster)
            end
        end
    end

    # 2b) adding nodes. For now, any adjacent nodes will do.
    # at some point, I would like to be able     


    i_cluster = 1
    subgraphs = []
    nodes_max_per_clusters =  nv(v_network) + nv(v_network) * (nv(s_network) - nv(v_network)) / nv(s_network)


    for cluster in clusters

        all_neighbors = Dict()
        for s_node in cluster
            for neigh in neighbors(s_network, s_node)
                if neigh ∉ cluster
                    if neigh ∉ keys(all_neighbors)
                        all_neighbors[neigh] = 1
                    else
                        all_neighbors[neigh] += 1
                    end
                end
            end
        end

        added = []
        while length(cluster) < nodes_max_per_clusters

            # ranking the neighbors²
            ranking = sort(collect(keys(all_neighbors)), by = x->all_neighbors[x], rev=true)
            # add the most connected neighbor
            push!(cluster, ranking[1])
            push!(added, ranking[1])
            delete!(all_neighbors, ranking[1])

            for neigh in neighbors(s_network, ranking[1])
                if neigh ∉ cluster
                    if neigh ∉ keys(all_neighbors)
                        all_neighbors[neigh] = 1
                    else
                        all_neighbors[neigh] += 1
                    end
                end
            end
        end
        sub_s_network = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")
        i_cluster += 1
        push!(subgraphs, sub_s_network)

        #write_added_nodes(s_network.graph, cluster, added, sub_s_network[][:name])
    end
        


    time_pretreatment = time() - time_beginning
    println("Well it took only $time_pretreatment to do all partitionning n all...")


    # 3 : Résoudre
    # We start by the most connected subgraph, because we are likely to spend more time on it.
    # It would be nice to take into account costs, capacities etc... But it looks a bit hard...


    
    ranking = sort(subgraphs, by = x->global_clustering_coefficient(x), rev=true)

    solving_time = 0
    i_cluster = 1
    time_beginning_solving = time()
    while solving_time < time_limit && i_cluster < length(subgraphs)+1
        sub_s_network = ranking[i_cluster]
        subinstance = Instance_Undir_VNE_1s(v_network, sub_s_network)

        current_model = Model(CPLEX.Optimizer)
        set_up_problem(subinstance, current_model)
        set_time_limit_sec(current_model, 1 )
        set_silent(current_model)

        optimize!(current_model)

        status_sol = primal_status(current_model)
    
        if status_sol == MOI.FEASIBLE_POINT
            println("New solution : $(objective_value(current_model))")
        else
            println("No solution found this time :'(")
        end

        i_cluster += 1

    end

    #visu_partitioning(s_network, partition)
    println("Found all em solution in just $(time() - time_beginning_solving)")

end



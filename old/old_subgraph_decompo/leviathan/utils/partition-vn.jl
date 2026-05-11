
includet("../../../utils/import_utils.jl")
includet("../../../utils/kahip_wrapper.jl")



function partition_vn(instance, nb_clusters)

    graph = instance.v_network.graph

    # 1 : Partitionner
    inbalance = 0.10
    println("$nb_clusters clusters to do, with inbalance $inbalance...")
    partition = partition_kahip(graph, nb_clusters, inbalance)
    clusters = [Vector{Int64}() for i in 1:nb_clusters]
    for s_node in vertices(graph)
        push!(clusters[partition[s_node]], s_node)
    end

    # 2 : Corriger

    for cluster in clusters
        simple_subgraph, vmap = induced_subgraph(graph, cluster)
        if !is_connected(simple_subgraph)
            components = connected_components(simple_subgraph)
            component_sorted = sort(components, by=x->length(x), rev=true)
            for subcluster in component_sorted[2:length(component_sorted)]
                #Let's add all those nodes to a (most) connected subgraph
                nodes_original = [vmap[node] for node in subcluster]
                subgraph_neighbors = zeros(Int, nb_clusters)
                for node in nodes_original
                    for neighbor in neighbors(graph, node)
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


    return clusters



end



function full_automatic_partition_vn(instance)
        
    v_network = instance.v_network 
    nb_max_part = ceil(Int, nv(v_network)/5)
    nb_min_part = max(2, floor(Int, nv(v_network)/15))
    println("Mmmmh let's do between $nb_min_part and $nb_max_part okay ?")


    best_score = 10000
    best_part = []
    # things to look at ? nb of connecting edges, average cluster coefficient of subgraphs...

    for nb_part in nb_min_part:nb_max_part

        println("let's do it for $nb_part !")

        graph = instance.v_network.graph

        # 1 : Partitionner
        inbalance = 0.10
        partition, edgecut = partition_kahip_edgecut(graph, nb_part, inbalance)
        clusters = [Vector{Int64}() for i in 1:nb_part]
        for s_node in vertices(graph)
            push!(clusters[partition[s_node]], s_node)
        end
    
        # 2 : Corriger
        for cluster in clusters
            simple_subgraph, vmap = induced_subgraph(graph, cluster)
            if !is_connected(simple_subgraph)
                components = connected_components(simple_subgraph)
                component_sorted = sort(components, by=x->length(x), rev=true)
                for subcluster in component_sorted[2:length(component_sorted)]
                    #Let's add all those nodes to a (most) connected subgraph
                    nodes_original = [vmap[node] for node in subcluster]
                    subgraph_neighbors = zeros(Int, nb_part)
                    for node in nodes_original
                        for neighbor in neighbors(graph, node)
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


        # 3 calcul score
        mean_cluster_factor = 0.1*nb_part
        for cluster in clusters
            simple_subgraph, vmap = induced_subgraph(graph, cluster)
            #mean_cluster_factor += ne(simple_subgraph) / (nv(simple_subgraph)^(3/2))
            # we compare the nb of edges to the number of edges in a perfect clique of that size uwu
            # that does penalize a bit larger graphs, 
            println("The global cluster coeff is : $(global_clustering_coefficient(simple_subgraph))")
            mean_cluster_factor += global_clustering_coefficient(simple_subgraph)
        end
        mean_cluster_factor = mean_cluster_factor/nb_part

        println("So for $nb_part")
        println("So we got an edge cut of $edgecut")
        println("And a mean cluster factor of $mean_cluster_factor")

        score = edgecut / (mean_cluster_factor)

        println("We got the score of $score")

        if score < best_score
            best_score = score
            best_part = clusters
        end

        
    end

    println("\n\nAt the end the one who won is $(length(best_part)) with $best_score")

end


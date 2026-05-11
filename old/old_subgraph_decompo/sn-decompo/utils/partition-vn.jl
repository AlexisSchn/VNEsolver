
includet("../../../utils/import_utils.jl")
includet("../../../utils/kahip_wrapper.jl")
includet("../../../utils/metis_wrapper.jl")


function partition_vn(instance, nb_clusters)

    graph = instance.v_network.graph
    
    # 1 : Partitionner
    inbalance = 0.1
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
                cluster_to_put_nodes_in = most_connected_subgraph[1]
                append!(clusters[cluster_to_put_nodes_in], nodes_original)
                for node in nodes_original
                    partition[node] = cluster_to_put_nodes_in
                end
                filter!(e->e∉nodes_original, cluster)
            end
        end
    end
    
    
    # equilibrating virtual subgraphs, to get all of them pretty close...
    target_size = nv(instance.v_network) / nb_clusters
    max_size = maximum(length(cluster) for cluster in clusters)
    min_size = minimum(length(cluster) for cluster in clusters)
    iter = 1
    while max_size > 1.2 * target_size || min_size < 0.8 * target_size

        if max_size > 1.2 * target_size # let's remove a node from the biggest subgraph
            oversized_clusters = findall(c -> length(c) > 1.2 * target_size, clusters)
            cluster_to_reduce = clusters[oversized_clusters[1]]

            neighboring_clusters = Dict()
            for node in cluster_to_reduce
                for neighbor in neighbors(graph, node)
                    if neighbor ∉ cluster_to_reduce
                        if partition[neighbor] ∈ keys(neighboring_clusters)
                            push!(neighboring_clusters[partition[neighbor]], node)
                        else
                            neighboring_clusters[partition[neighbor]] = [node]
                        end
                    end
                end
            end
            #println("hum sooo... $neighboring_clusters")
            smallest_cluster_around = argmin(i -> length(clusters[i]), keys(neighboring_clusters))
            #println("Well, let's add it to $smallest_cluster_around !")
            node_to_change_cluster = neighboring_clusters[smallest_cluster_around][1]
            deleteat!(cluster_to_reduce, findfirst(==(node_to_change_cluster), cluster_to_reduce)) 
            push!(clusters[smallest_cluster_around], node_to_change_cluster)
            partition[node_to_change_cluster] = smallest_cluster_around
            # put a node in the neighbor cluster with fewer nodes



        else # let's add a node to the smallest subgraph


            # do stuff..


        end




        # Update max size
        max_size = maximum(length(cluster) for cluster in clusters)
        min_size = minimum(length(cluster) for cluster in clusters)

        if iter == 1
            println("Have to do extra balancing of the subgraphs, due to KaHIP...")
        end
        iter += 1
        if iter > 30
            println("Coudn't balance the subgraphs well... Sorry.")
            break
        end
    end


    return clusters

end




function partition_vn_metis(instance, nb_clusters)

    graph = instance.v_network.graph
    
    # 1 : Partitionner
    println("$nb_clusters clusters to do... METIS TIME!")
    partition = partition_metis(graph, nb_clusters, 1.)
    clusters = [Vector{Int64}() for i in 1:nb_clusters]
    for s_node in vertices(graph)
        push!(clusters[partition[s_node]], s_node)
    end

    # 2 : Corriger
    for cluster in clusters
        simple_subgraph, vmap = induced_subgraph(graph, cluster)
        if !is_connected(simple_subgraph)
            println("PUTAIN C PAS CONNECTE")
        end
    end

    # To do : So post-treatment, because, it's still NOT BALANCED!!!

    return clusters
end

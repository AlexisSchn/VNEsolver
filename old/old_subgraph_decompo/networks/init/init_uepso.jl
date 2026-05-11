
includet("../../../utils/import_utils.jl")
includet("../../../utils/metis_wrapper.jl")

includet("../../../heuristics/uepso.jl")
#includet("../../../heuristics/pso2.jl")
includet("../../../heuristics/pso.jl")
includet("../../../heuristics/mepso.jl")

includet("../utils/utils-subgraphdecompo.jl")


# Nb columns is number of cols PER virtual subgraph !
function init_uepso(instance, vn_decompo, nb_columns)


    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    # Some way to partition the sn network into nb subgraph parts
    s_subgraphs = get_sn_decompo_kahip(s_network, nb_columns, 20)

    

    # Get the mappings (columns)
    mappings_result = Dict()
    
    for v_subgraph in vn_decompo.subgraphs

        mappings = []

        for s_subgraph in s_subgraphs

            sub_instance = Instance(v_subgraph.graph, s_subgraph.graph)
            
            sub_mapping, cost = solve_mepso(sub_instance; nb_particle=30, nb_iter=50, time_max=0.25, print_things=false)

            if isnothing(sub_mapping)
                continue
            end
            true_cost = 0
            node_placement = []
            for v_node in vertices(v_subgraph.graph)
                real_s_node = s_subgraph.nodes_of_main_graph[sub_mapping.node_placement][v_node]
                append!(node_placement, real_s_node)
                true_cost += s_network[real_s_node][:cost]
            end
    
    
            edge_routing = Dict()
            for v_edge in edges(v_subgraph.graph)
                used_edges = []
                for s_edge in sub_mapping.edge_routing[v_edge].edges
                    real_s_edge = get_edge(s_network_dir, s_subgraph.nodes_of_main_graph[src(s_edge)], s_subgraph.nodes_of_main_graph[dst(s_edge)])
                    push!(used_edges, real_s_edge)
                    true_cost += s_network_dir[src(real_s_edge), dst(real_s_edge)][:cost]
                end
                edge_routing[v_edge] = order_path(s_network_dir, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
            end
            real_sub_mapping = Mapping(v_subgraph.graph, s_network_dir, node_placement, edge_routing)
    
            push!(mappings, real_sub_mapping)

        end

        mappings_result[v_subgraph] = mappings
    end


    return mappings_result

end






# computing the sn decompo. It depends on the subvirtual network...
function get_sn_decompo_kahip(s_network, nb_clusters, nb_nodes_per_clusters)


    # 1 : Partitionner
    inbalance = 0.10
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

    i_cluster = 1
    subgraphs = []

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

        nb_nodes_with_capacity=0
        for s_node in cluster
            if s_network[s_node][:cap] >= 1
                nb_nodes_with_capacity+=1
            end
        end

        while nb_nodes_with_capacity < nb_nodes_per_clusters

            # ranking the neighbors
            ranking = sort(collect(keys(all_neighbors)), by = x->all_neighbors[x], rev=true)
            # add the most connected neighbor
            new_s_node = ranking[1]
            if s_network[new_s_node][:cap] >= 1
                nb_nodes_with_capacity+=1
            end
            push!(cluster, new_s_node)
            push!(added, new_s_node)
            delete!(all_neighbors, new_s_node)

            for neigh in neighbors(s_network, new_s_node)
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
        s_subg = Subgraph(sub_s_network, cluster)
        i_cluster += 1
        push!(subgraphs, s_subg)

    end



    return subgraphs

end









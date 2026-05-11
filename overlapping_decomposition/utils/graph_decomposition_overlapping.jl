using Graphs, MetaGraphsNext


### === Structs
struct NetworkDecompositionOverlapping
    subgraphs
    v_nodes_assignment
    v_edges_master
    overlapping_nodes
end

struct Subgraph
    graph
    nodes_of_main_graph
end


function print_stuff_subgraphs(original_graph, subgraphs)

    println("For $original_graph, there is $(length(subgraphs)) subgraphs:")
    for subgraph in subgraphs
        println("       $(subgraph.graph[][:name]) with $(nv(subgraph.graph)) nodes and $(ne(subgraph.graph)) edges")
    end
    
end


function set_up_decompo_overlapping(instance, node_partitionning)

    vn = instance.v_network

    # Here, node assigment is a dict, I don't know why. And node assigment ofa  node is alsoa dic,
    # because it can be in several sub-networks.
    node_assignment = Dict()
    for v_node in vertices(vn)
        node_assignment[v_node] = Dict()
    end

    # getting the subgraphs and the node assignment
    subgraphs = []
    for (i_subgraph, v_nodes) in enumerate(node_partitionning)
        subgraph = Subgraph(my_induced_subgraph(vn, v_nodes, "subgraph_$i_subgraph"), v_nodes)
        for (i_node, v_node) in enumerate(v_nodes)
            node_assignment[v_node][subgraph] = i_node
        end
        push!(subgraphs, subgraph)
    end


    # finding and removing overlapping edges!
    for v_edge in edges(vn)
        common_subgraph = collect(keys(node_assignment[src(v_edge)]) ∩ keys(node_assignment[dst(v_edge)]))
        if length(common_subgraph) ≥ 2
            println("V edge $v_edge is in more than one subgraph: $common_subgraph")
            size_subg = [ne(subgraph.graph) for subgraph in common_subgraph]
            idx_ranking = sortperm(size_subg)
            for idx_subgraph in idx_ranking[2:end]
                subgraph = common_subgraph[idx_subgraph]
                src_in_subgraph = node_assignment[src(v_edge)][subgraph]
                dst_in_subgraph = node_assignment[dst(v_edge)][subgraph]
                rem_edge!(subgraph.graph, src_in_subgraph, dst_in_subgraph)
                if !is_connected(subgraph.graph)
                    println("Well, the graph is not connected anymore... That sucks... Fix this or include overlapping edges...")
                    # remove one of the two nodes from  the subgraph ?
                    if degree(subgraph.graph, src_in_subgraph) == 0
                        rem_vertex!(subgraph.graph, src_in_subgraph)
                        filter!(x -> x != src(v_edge), subgraph.nodes_of_main_graph)
                        delete!(node_assignment[src(v_edge)], subgraph)
                        println("I removed $(src(v_edge)) from $subgraph !")
                    end
                    if degree(subgraph.graph, dst_in_subgraph) == 0
                        rem_vertex!(subgraph.graph, dst_in_subgraph)
                        filter!(x -> x != dst(v_edge), subgraph.nodes_of_main_graph)
                        delete!(node_assignment[dst(v_edge)], subgraph)
                        println("I removed $(dst(v_edge)) from $subgraph !")
                    end
                    if !is_connected(subgraph.graph)
                        println("WTF it's still not connected??")
                    end
                end
            end
            println("I removed the edge: $common_subgraph")
        end
    end


    # finding overlapping nodes
    v_node_overlapping = Dict()
    for v_node in vertices(vn)
        if length(keys(node_assignment[v_node])) > 1
            v_node_overlapping[v_node] = keys(node_assignment[v_node])
            println("$v_node is overlapping !")
        end
    end

    
    # finding out the master virtual edges
    v_edge_master = [] 
    for v_edge in edges(vn)
        in_master = true
        for subgraph_src in keys(node_assignment[src(v_edge)])
            for subgraph_dst in keys(node_assignment[dst(v_edge)])
                if subgraph_src == subgraph_dst
                    in_master = false
                end
            end
        end
        if in_master
            push!(v_edge_master, v_edge)
        end
    end

    vn_decompo = NetworkDecompositionOverlapping(subgraphs, node_assignment, v_edge_master, v_node_overlapping)

    return vn_decompo
end



# if you wanna also tell which edges to take
# sub_vns is a list of dict, dict includes "nodes" and "edges"
function set_up_decompo_overlapping_more_info(instance, sub_vns)

    vn = instance.v_network

    # Here, node assigment is a dict, I don't know why. And node assigment of a node is also a dict,
    # because it can be in several sub-networks.
    node_assignment = Dict()
    for v_node in vertices(vn)
        node_assignment[v_node] = Dict()
    end

    # getting the subgraphs and the node assignment
    edges_in_subgraphs = []
    subgraphs = []
    for (i_subgraph, sub_vn) in enumerate(sub_vns)
        println("Subgrpah $i_subgraph: $sub_vn")
        v_nodes = sub_vn["nodes"]
        subgraph = Subgraph(my_induced_subgraph(vn, v_nodes, "subgraph_$i_subgraph"), v_nodes)
        for (i_node, v_node) in enumerate(v_nodes)
            node_assignment[v_node][subgraph] = i_node
        end
        for v_edge in edges(subgraph.graph)
            remove_it = true
            for v_edge_user in sub_vn["edges"]
                if (src(v_edge_user) == v_nodes[src(v_edge)]) && (dst(v_edge_user) == v_nodes[dst(v_edge)])
                    remove_it = false
                end
                if (dst(v_edge_user) == v_nodes[src(v_edge)]) && (src(v_edge_user) == v_nodes[dst(v_edge)])
                    remove_it = false
                end 
            end
            if remove_it
                rem_edge!(subgraph.graph, src(v_edge), dst(v_edge))
            end
        end
        push!(subgraphs, subgraph)
        append!(edges_in_subgraphs, sub_vn["edges"])
    end

    # Maybe put something to check for overlapping edges?



    # finding overlapping nodes
    v_node_overlapping = Dict()
    for v_node in vertices(vn)
        if length(keys(node_assignment[v_node])) > 1
            v_node_overlapping[v_node] = keys(node_assignment[v_node])
            println("$v_node is overlapping !")
        end
    end

    
    # finding out the master virtual edges
    # a bit more annoying here: must be sure 
    v_edge_master = [] 
    for v_edge in edges(vn)
        in_master = true
        for v_edge_subg in edges_in_subgraphs
            if (src(v_edge) == src(v_edge_subg)) && (dst(v_edge) == dst(v_edge_subg))
                in_master = false
            end
            if (dst(v_edge) == src(v_edge_subg)) && (src(v_edge) == dst(v_edge_subg))
                in_master = false
            end
        end
        if in_master
            push!(v_edge_master, v_edge)
        end
    end

    vn_decompo = NetworkDecompositionOverlapping(subgraphs, node_assignment, v_edge_master, v_node_overlapping)

    return vn_decompo

end

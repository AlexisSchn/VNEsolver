using Graphs, MetaGraphsNext

includet("../utils/import_utils.jl")
includet("shortest-path-routing.jl")
includet("optimal_routing.jl")


# Implementation of : Virtual Network Embedding Based on the Degree and Clustering Coefficient Information
# I did some slight modifications. It would be better to get a clean version...

function solve_two_stage_DCC(instance)


    v_network = instance.v_network
    s_network = instance.s_network
    s_network = instance.s_network



    # --- STEP 1 --- RANK THE NODES

    # Virtual node ranking: the first one is the most central. Then, it is the closest the it, in term of shortest path.
    v_nodes_centrality = closeness_centrality(v_network)
    root = findmax(v_nodes_centrality)[2]
    distances = desopo_pape_shortest_paths(instance.v_network, root).dists 
    v_nodes_scores = [ v_nodes_centrality[v_node] - distances[v_node] * 10 for v_node in vertices(v_network)]
    v_nodes_ordered = sortperm(v_nodes_scores; rev=true)
    
    # Substrate node ranking
    s_nodes_capacities_around = [s_network[s_node][:cap] * sum(s_network[src(s_edge), dst(s_edge)][:cap] 
                        for s_edge in get_out_edges(s_network, s_node) ∪ get_in_edges(s_network, s_node))
                            for s_node in vertices(s_network)
    ]
    
    s_node_centrality = closeness_centrality(s_network)

    s_nodes_scores = [s_node_centrality[s_node] * s_nodes_capacities_around[s_node] for s_node in vertices(s_network)]
    s_nodes_ordered = sortperm(s_nodes_scores; rev=true)




    # --- STEP 2 --- PLACE THE VIRTUAL NODES
    v_node_placement = Dict()
    v_node_placed = []
    s_node_used = []
    v_node_placement_cost = 0

    # init: we place the most central v nodes on the most central s nodes
    v_node_placement[v_nodes_ordered[1]] = s_nodes_ordered[1]
    push!(v_node_placed, v_nodes_ordered[1])
    push!(s_node_used, s_nodes_ordered[1])
    v_node_placement_cost += get_attribute_node(s_network, s_nodes_ordered[1], :cost)

    for v_node in v_nodes_ordered[2:end]
        # find the sum of the distance between a s_node and the v_node_placements of neighbors of v_node that have been placed. 
        # here is the improvement, in the article they only do one neighbor, based on the tree
        s_nodes_of_neighbors = []
        for neighbor in neighbors(v_network, v_node)
            if neighbor ∈ v_node_placed
                push!(s_nodes_of_neighbors, v_node_placement[neighbor])
            end
        end

        current_distances_s_nodes = [sum(desopo_pape_shortest_paths(instance.s_network, placement_neighbor).dists[s_node] for placement_neighbor in s_nodes_of_neighbors) for s_node in vertices(s_network)]
        current_s_node_scores = [s_nodes_scores[s_node] - current_distances_s_nodes[s_node] * 10 for s_node in vertices(s_network)]

        current_s_nodes_ranked = sortperm(current_s_node_scores; rev=true)
        keep_on = true
        idx_s_node = 1
        while keep_on
            s_node = current_s_nodes_ranked[idx_s_node]
            if s_node ∉ s_node_used && s_network[s_node][:cap] >= 1 
                v_node_placement[v_node] = s_node
                push!(v_node_placed, v_node)
                push!(s_node_used, s_node)
                keep_on = false
                v_node_placement_cost += get_attribute_node(s_network, s_node, :cost)
            end
            idx_s_node += 1
        end
    end
    #println("v_node_placement of nodes: $v_node_placement")



    # STEP 3 --- ROUTE THE VIRTUAL EDGES
    v_edge_routing, routing_cost = shortest_path_routing(instance, v_node_placement)
    #v_edge_routing, routing_cost = optimal_routing(instance, v_node_placement)


    println("Overall, there is a cost of $routing_cost for routing, $v_node_placement_cost for placement")

    return (v_node_placement_cost + routing_cost)
end






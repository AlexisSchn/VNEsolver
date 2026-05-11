using Graphs, MetaGraphsNext
includet("../utils/import_utils.jl")



# from an edge placement, we compute the edge routing, by looking at shortest paths.
# this is a common algorithm from the literature.
# usually, they also do some type of load balancing - but since we are doing offline vne, we are not interested in that here.
function shortest_path_routing(instance, v_node_placement)


    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir
    s_network_dir_copy = deepcopy(instance.s_network_dir)

    edge_routing = Dict()

    # dstmax matrix
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(instance.s_network)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(instance.s_network, s_edge, :cost)
        distmx[dst(s_edge), src(s_edge)] = get_attribute_edge(instance.s_network, s_edge, :cost)
    end

    overall_cost = 0
    for v_edge in edges(v_network)
        s_src = v_node_placement[src(v_edge)]
        s_dst = v_node_placement[dst(v_edge)]
        shortest_path = a_star(s_network_dir_copy, s_src, s_dst, distmx)


        if shortest_path == []
            #println("No shortest path found: the graph is full!")
            #println("I had the following routing: $edge_routing")
            return Dict(), 99999999
        end

        edge_routing[v_edge] = order_path(s_network_dir, shortest_path, s_src, s_dst) 

        for s_edge in shortest_path
            set_attribute_edge(s_network_dir_copy, s_edge, :cap,  get_attribute_edge(s_network_dir_copy, s_edge, :cap)-1)
            set_attribute_edge(s_network_dir_copy, get_reverse_edge(s_network_dir_copy, s_edge), :cap,  get_attribute_edge(s_network_dir_copy, get_reverse_edge(s_network_dir_copy, s_edge), :cap)-1)
            if get_attribute_edge(s_network_dir_copy, s_edge, :cap) <= 0
                #println("Well it is time to stop using $s_edge kids")
                distmx[src(s_edge), dst(s_edge)] = 0
                distmx[dst(s_edge), src(s_edge)] = 0
                rem_edge!(s_network_dir_copy, src(s_edge), dst(s_edge))
                rem_edge!(s_network_dir_copy, dst(s_edge), src(s_edge))
            end
        end

        overall_cost += edge_routing[v_edge].cost
    end

    #println("Well, here is a nice edge routing: $edge_routing of cost $overall_cost")
    
    return edge_routing, overall_cost
end


# from an edge placement, we compute the edge routing, by looking at shortest paths.
# Here, we allow some violation of the capacities. Indeed, otherwise, the heuristic algorithms that relies on this can not work at all.
#! Looks like it doesnt help at all.
function shortest_path_routing_nicer(instance, v_node_placement)


    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir
    s_network_dir_copy = deepcopy(instance.s_network_dir)

    edge_routing = Dict()

    # dstmax matrix
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(instance.s_network)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(instance.s_network, s_edge, :cost)
        distmx[dst(s_edge), src(s_edge)] = get_attribute_edge(instance.s_network, s_edge, :cost)
    end

    overall_cost = 0
    for v_edge in edges(v_network)
        s_src = v_node_placement[src(v_edge)]
        s_dst = v_node_placement[dst(v_edge)]
        shortest_path = a_star(s_network_dir_copy, s_src, s_dst, distmx)

        edge_routing[v_edge] = order_path(s_network_dir_copy, shortest_path, s_src, s_dst) 

        for s_edge in shortest_path  # Here, is there is not enough capacity, we "allow" the violation (ie. we don't remove the edge), but the price of the edge will be very high
            if get_attribute_edge(s_network_dir_copy, s_edge, :cap) > 1
                set_attribute_edge(s_network_dir_copy, s_edge, :cap,  get_attribute_edge(s_network_dir_copy, s_edge, :cap)-1)
                set_attribute_edge(s_network_dir_copy, get_reverse_edge(s_network_dir_copy, s_edge), :cap,  get_attribute_edge(s_network_dir_copy, get_reverse_edge(s_network_dir_copy, s_edge), :cap)-1)
            else

                if get_attribute_edge(s_network_dir_copy, s_edge, :cost) < 500
                    set_attribute_edge(s_network_dir_copy, s_edge, :cost,  500)
                    set_attribute_edge(s_network_dir_copy, get_reverse_edge(s_network_dir_copy, s_edge), :cost, 500)
                end
            end
        end

        overall_cost += edge_routing[v_edge].cost
    end

    #println("Well, here is a nice edge routing: $edge_routing of cost $overall_cost")
    
    return edge_routing, overall_cost
end

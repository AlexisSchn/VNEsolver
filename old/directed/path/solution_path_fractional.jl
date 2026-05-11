

struct Mapping_Path_Fractional
    virtual_network
    substrate_network
    node_placement
    node_placement_cost
    edge_routing
    edge_routing_cost
end


function Mapping_Path_Fractional(v_network, s_network, node_placement, edge_routing)
    node_placement_cost = 0
    for v_node in vertices(v_network)
        for s_node in vertices(s_network)
            node_placement_cost += node_placement[v_node][s_node] * s_network[s_node][:cost] * v_network[v_node][:dem]
        end
    end

    edge_routing_cost = 0
    for v_edge in edges(v_network)
        for s_path in keys(edge_routing[v_edge])
            edge_routing_cost += s_path.cost * edge_routing[v_edge][s_path] * v_network[src(v_edge), dst(v_edge)][:dem]
        end
    end

    return Mapping_Path_Fractional(v_network, s_network, node_placement, node_placement_cost, edge_routing, edge_routing_cost)
end


function Base.show(io::IO, mapping::Mapping_Path_Fractional)
    println(io, "Node placement of cost : " * string(round(mapping.node_placement_cost * 1000) / 1000))
    for v_node in vertices(mapping.virtual_network)
        println(io, string(v_node) * " : ")
        for s_node in vertices(mapping.substrate_network)
            if mapping.node_placement[v_node][s_node] > 0.001
                println(io, "       " * string(s_node) * " : " * string(round(mapping.node_placement[v_node][s_node] * 1000) / 1000))
            end
        end
    end

    println(io, "Edge routing of cost: " * string(round(mapping.edge_routing_cost * 1000) / 1000))
    for v_edge in edges(mapping.virtual_network)
        println(io, string(v_edge) * " : ")
        for s_path in keys(mapping.edge_routing[v_edge])
            if mapping.edge_routing[v_edge][s_path] > 0.001
                println(io, "       "  * string(round(mapping.edge_routing[v_edge][s_path]*1000)/1000) * " : " * string(s_path))
            end
        end
    end
end


function print_solution_path_fractional(mappings)
    println("Printing fractional solution of path solving:\n")
    i_mapping = 1
    for mapping in mappings
        println("Printing mapping of virtual graph " * string(i_mapping))
        println(mapping)
        i_mapping += 1
    end
    print("\n")
end


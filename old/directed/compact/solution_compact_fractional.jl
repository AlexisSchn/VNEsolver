### MAPPING ###

struct MappingCompactFractional
    virtual_network
    substrate_network
    node_placement
    node_placement_cost
    edge_routing
    edge_routing_cost
end

function MappingCompactFractional(v_network, s_network, node_placement, edge_routing)
    node_placement_cost = 0
    for v_node in vertices(v_network)
        for s_node in vertices(s_network)
            node_placement_cost += node_placement[v_node][s_node] * s_network[s_node][:cost] * v_network[v_node][:dem]
        end
    end

    edge_routing_cost = 0
    for v_edge in edges(v_network)
        for s_edge in edges(s_network)
            edge_routing_cost += edge_routing[v_edge][s_edge] * s_network[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem]
        end
    end

    return MappingCompactFractional(v_network, s_network, node_placement, node_placement_cost, edge_routing, edge_routing_cost)
end



function Base.show(io::IO, mapping::MappingCompactFractional)
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
        for s_edge in edges(mapping.substrate_network)
            if mapping.edge_routing[v_edge][s_edge] > 0.001
                println(io, "       "  * string(round(mapping.edge_routing[v_edge][s_edge]*1000)/1000) * " : " * string(s_edge))
            end
        end
    end
end


function print_solution_compact_fractional(mappings)
    println("Printing compact fractional solution :\n")
    i_mapping = 1
    for mapping in mappings
        println("Printing mapping of virtual graph " * string(i_mapping))
        println(mapping)
        i_mapping += 1
    end
    print("\n")
end


function solution_checker(instance, mappings)
    solution_is_valid = true

    # Capacité noeuds
    for s_node in vertices(instance.s_network)
        dem_totale = 0
        i_mapping = 1
        for v_network in instance.v_networks
            current_mapping = mappings[i_mapping]
            for v_node in vertices(v_network)
                if current_mapping.node_placement[v_node] == s_node
                    dem_totale += v_network[v_node][:dem]
                end
            end
            i_mapping += 1
        end
        if dem_totale > instance.s_network[s_node][:cap]
            println("Capacity of node " * string(s_node) * " not respected : " * string(dem_totale) * " used")
            solution_is_valid = false
        end
    end
    
    
    # Capacité arretes
    for s_edge in edges(instance.s_network)
        dem_totale = 0
        i_mapping = 1
        for v_network in instance.v_networks
            current_mapping = mappings[i_mapping]
            for v_edge in edges(v_network)
                if s_edge in current_mapping.edge_routing[v_edge].edges
                    dem_totale += v_network[src(v_edge), dst(v_edge)][:dem]
                end
            end
            i_mapping += 1
        end
        cap_edge = instance.s_network[src(s_edge), dst(s_edge)][:cap]
        if dem_totale > cap_edge
            println("Capacity of edge " * string(s_edge) * " not respected : " * string(dem_totale) * " used for " * string(cap_edge))
            solution_is_valid = false
        end
    end
    
    # Path departure and arrival at right node ? todo
    


    if solution_is_valid
        println("Solution is valid")
    else
        throw("Error : invalid solution")
    end
end


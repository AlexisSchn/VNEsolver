### MAPPING ###

struct Mapping
    virtual_network
    substrate_network
    node_placement
    edge_routing
end





function get_cost_placement(mapping)
    v_network = mapping.virtual_network
    s_network = mapping.substrate_network
    node_placement_cost = 0
    for v_node in vertices(v_network)
        node_placement_cost += s_network[mapping.node_placement[v_node]][:cost]
    end

    return node_placement_cost
end


function get_usage_placement(mapping)

    s_node_usage = []
    for s_node in vertices(s_network)
        usage = 0
        for v_node in vertices(v_network)
            if node_placement[v_node] == s_node
                usage += 1
            end
        end
        push!(s_node_usage, usage)
    end

    return s_node_usage
end


function get_cost_routing(mapping)
    v_network = mapping.virtual_network

    edge_routing_cost = 0
    for v_edge in edges(v_network)
        path = mapping.edge_routing[v_edge]
        edge_routing_cost += path.cost
    end
    return edge_routing_cost
end

function get_usage_routing(mapping)

    s_edge_usage = Dict()
    for s_edge in edges(s_network)
        usage = 0
        for v_edge in edges(v_network)
            if s_edge in edge_routing[v_edge].edges
                usage += v_network[src(v_edge), dst(v_edge)][:dem]
            end
        end
        s_edge_usage[s_edge] = usage
    end

    return s_edge_usage 
end




function Base.show(io::IO, mapping::Mapping)
    print(io, "Mapping of $(mapping.virtual_network) on $(mapping.substrate_network) with placement $(mapping.node_placement)")
end



function solution_checker(instance, mapping)
    solution_is_valid = true

    # Capacité noeuds
    for s_node in vertices(instance.substrate_network)
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
    for s_edge in edges(instance.substrate_network)
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
        cap_edge = instance.s_network[src(s_edge), dst(d_edge)][:cap]
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

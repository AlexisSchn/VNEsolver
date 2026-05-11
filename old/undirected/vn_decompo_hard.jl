using Revise

using Graphs, MetaGraphsNext

includet("../../utils/import_utils.jl")
includet("../directed/vn_decompo_full/VNDecompoFull.jl")

using .NetworkDecompositionFull

function solve_undir_vndecompo_colge_hard(instance, node_partitionning)

    # Construct directed instance
    instance_dir = get_directed_instance(instance)
    println("Instance directed constructed, starting colge...")
    
    solve_dir_vn_decompo(instance_dir, node_partitionning);

    

end



function get_directed_instance(instance)
    s_network_dir = get_directed_sn(instance.s_network)
    v_networks_dir = []
    for v_network in instance.v_networks
        push!(v_networks_dir, get_directed_vn(v_network))
    end
    return InstanceVNE(v_networks_dir, s_network_dir)
end


function get_directed_sn(s_network)
    sn_dir = MetaGraph(
        DiGraph(),
        Int,
        Dict,
        Dict,
        Dict(:name => s_network[][:name], :type => s_network[][:type], :directed => true)
    )
    
    for s_node in vertices(s_network)
        add_vertex!(sn_dir, s_node, s_network[s_node])
    end
    # Dividing by 2 the edge cost, because the path will need to be done both ways !
    for s_edge in edges(s_network)
        add_edge!(sn_dir, src(s_edge), dst(s_edge), copy(s_network[src(s_edge), dst(s_edge)]))
        sn_dir[src(s_edge), dst(s_edge)][:cost] = s_network[src(s_edge), dst(s_edge)][:cost] / 2
        add_edge!(sn_dir, dst(s_edge), src(s_edge), copy(s_network[src(s_edge), dst(s_edge)]))
        sn_dir[dst(s_edge), src(s_edge)][:cost] = s_network[src(s_edge), dst(s_edge)][:cost] / 2
    end
    return sn_dir
end


function get_directed_vn(v_network)
    vn_dir = MetaGraph(
        DiGraph(),
        Int,
        Dict,
        Dict,
        Dict(:name => v_network[][:name], :type => v_network[][:type], :directed => true)
    )
    
    for v_node in vertices(v_network)
        add_vertex!(vn_dir, v_node, v_network[v_node])
    end
    
    for v_edge in edges(v_network)
        add_edge!(vn_dir, src(v_edge), dst(v_edge), copy(v_network[src(v_edge), dst(v_edge)]))
        add_edge!(vn_dir, dst(v_edge), src(v_edge), copy(v_network[src(v_edge), dst(v_edge)]))
    end
    return vn_dir
end


function make_solution_undir(instance_undir, mappings_dir)

    mappings_undir = []
    s_network_undir = instance_undir.s_network

    for i_vn in 1:length(mappings_dir)
        v_network = instance_undir.v_networks[i_vn]
        mapping_vn_dir = mappings_dir[i_vn]
        node_placement = mapping_vn_dir.node_placement
        edge_routing = Dict()
        for v_edge in edges(v_network)
            path_dir = mapping_vn_dir.edge_routing[v_edge]
            if path_dir.src == path_dir.dst
                edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
            else
                edges_used_undir = []
                for edge_dir in path_dir.edges
                    push!(edges_used_undir, get_edge(s_network_undir, src(edge_dir), dst(edge_dir)))
                end
                edge_routing[v_edge] = order_path(s_network_undir, edges_used_undir, node_placement[src(v_edge)], node_placement[dst(v_edge)])
            end
        end
        push!(mappings_undir, Mapping(v_network, instance_undir.s_network, node_placement, edge_routing) )
    end


    return mappings_undir
end
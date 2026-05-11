includet("../../utils/import_utils.jl")
includet("../directed/compact/compact_formulation.jl")
includet("../directed/vn_decompo/vn_decompo.jl")
#includet("../directed/cuts/cuts.jl")





function solve_undir_compact_dir(instance, time_solver = 30, silent = false)
    # Construct directed instance
    instance_dir = get_directed_instance(instance)

    # solve 
    # problem = set_up_compact_model_furobi(instance_dir, true, true, true)
    problem = set_up_compact_model(instance_dir, true, true, true)

    #@constraint(problem.model, sum(problem.model[:y][v_network, v_edge, s_edge] for v_network in instance_dir.v_networks for v_edge in edges(v_network) for s_edge in edges(instance_dir.s_network)) >= 70)

    print("Starting solving... ")
    set_time_limit_sec(problem.model, time_solver)
    
    if silent
        set_silent(problem.model)
    end
    optimize!(problem.model)
    println("done. Solving state: " * string(termination_status(problem.model)) * ", obj value: " * 
            string(objective_value(problem.model)) * ", bound value: " * string(objective_bound(problem.model)))

    # Get the solution
    x_values = value.(problem.model[:x])
    y_values = value.(problem.model[:y])
    mappings_dir = get_solution(instance_dir, x_values, y_values)

    mappings_undir = make_solution_undir(instance, mappings_dir)

    #println("Mappings dir : " * string(mappings_dir))
    #println("Mappings undir : " * string(mappings_undir))

    # get correct solution

    return mappings_undir
end




function solve_undir_compact_fractional_dir(instance, time_solver = 30, silent = false)
    # Construct directed instance
    instance_dir = get_directed_instance(instance)

    # solve 
    # problem = set_up_compact_model_furobi(instance_dir, true, true, true)
    mappings = solve_directed_compact_fractional(instance_dir, true, true)
    for mapping in mappings
        println(mapping)
    end
    return mappings_undir
end



function solve_undir_vndecompo_colge(instance, node_partitionning)

    # Construct directed instance
    instance_dir = get_directed_instance(instance)
    println("Instance directed constructed, starting colge...")
    
    mappings = vn_decompo(instance_dir, node_partitionning)
    println("OKAY LET'S LOOK AT THE MAPPINGS: ")
    #println(mappings[1])
    

end



function solve_undir_cuts_cycles(instance, cycles)
    # Construct directed instance
    instance_dir = get_directed_instance(instance)
 
    solve_directed_cuts_cycles(instance_dir, cycles, true, true, 30, false)
end





function solve_undir_compact_subgraphcst(instance, node_partitionning, time_solving = 30)
    
    instance_dir = get_directed_instance(instance)

    mappings = solve_directed_compact_with_subgraphs_constraints(instance_dir, node_partitionning, time_solving)

    print(mappings)
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
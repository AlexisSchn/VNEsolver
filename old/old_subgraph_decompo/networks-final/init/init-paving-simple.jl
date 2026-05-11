
using Graphs, MetaGraphsNext
includet("../../../heuristics/mepso.jl")
includet("../../../compact/compact_plus.jl")
includet("../../../heuristics/local-search-heuristic.jl")


function find_submappings_simple(instance, vn_decompo, sn_subgraphs; solver="mepso", nb_columns=200)


    s_network = instance.s_network
    s_network_dir = instance.s_network_dir
    v_network = instance.v_network
    vn_subgraphs = vn_decompo.subgraphs 
    mappings = []
    mappings_per_subgraph = Dict()
    for v_subgraph in vn_subgraphs
        mappings_per_subgraph[v_subgraph] = []
    end
    nb_substrate_subgraphs = length(sn_subgraphs)

    


    # Base shortest paths
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network_dir)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir, s_edge, :cost)
    end
    base_shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)

        
    # All vn on all sn
    while length(mappings) < nb_columns

        # Associate subvn to a random subsn
        used_sub_s_network = []
        assignment_virtual_substrate_subgraphs = Dict()
        for v_subgraph in vn_subgraphs
            found = false
            while !found
                i_subgraph = rand(1:nb_substrate_subgraphs)
                if i_subgraph ∉ used_sub_s_network
                    push!(used_sub_s_network, i_subgraph)
                    found = true
                    assignment_virtual_substrate_subgraphs[v_subgraph] = sn_subgraphs[i_subgraph]
                end
            end
        end



        for v_subgraph in vn_subgraphs
            s_subgraph = assignment_virtual_substrate_subgraphs[v_subgraph]
            sub_instance= Instance(v_subgraph.graph, s_subgraph.graph)

            if solver=="mepso"
                sub_mapping, cost = solve_mepso(sub_instance; nb_particle=25, nb_iter=50, time_max=0.2, print_things=false)
            elseif solver=="local-search"
                result=solve_local_search(sub_instance; nb_particle=25, nb_local_search = 50)
                sub_mapping = result["mapping"]
            elseif solver=="milp"
                result=solve_compact_ffplus(sub_instance)
                sub_mapping = result["mapping"]
            else
                println("I don't know your solver. using MILP.")
                result=solve_compact_ffplus(sub_instance)
                sub_mapping = result["mapping"]
            end

            if isnothing(sub_mapping)
                continue
            end
            
            node_placement = []
            for v_node in vertices(v_subgraph.graph)
                original_s_node = s_subgraph.nodes_of_main_graph[sub_mapping.node_placement][v_node]
                append!(node_placement, original_s_node)
            end


            edge_routing = Dict()
            for v_edge in edges(v_subgraph.graph)
                used_edges = []
                for s_edge in sub_mapping.edge_routing[v_edge].edges
                    real_s_edge = get_edge(s_network_dir, s_subgraph.nodes_of_main_graph[src(s_edge)], s_subgraph.nodes_of_main_graph[dst(s_edge)])
                    push!(used_edges, real_s_edge)
                end
                edge_routing[v_edge] = order_path(s_network_dir, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
            end


            real_sub_mapping = Mapping(v_subgraph.graph, s_network_dir, node_placement, edge_routing)
    
            push!(mappings, real_sub_mapping)
            push!(mappings_per_subgraph[v_subgraph], real_sub_mapping)

        end
        
        print("We have $(length(mappings)) mappings! ")

    end

    println("We have obtained $(length(mappings)) mappings!")
    return mappings_per_subgraph    
end


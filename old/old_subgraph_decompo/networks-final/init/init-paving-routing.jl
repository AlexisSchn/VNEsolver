
using Graphs, MetaGraphsNext
includet("../heuristics/mepso-enhanced.jl")
includet("../heuristics/ilp-routing-cost.jl")
includet("../heuristics/local-search-routing.jl")


function find_submappings_routing(instance, vn_decompo, sn_subgraphs; solver="mepso", nb_columns=200)


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

        temporary_placement = zeros(Int, nv(v_network))
        
        for v_subgraph in vn_subgraphs
            s_subgraph = assignment_virtual_substrate_subgraphs[v_subgraph]
            cluster = s_subgraph.nodes_of_main_graph
            capacities_score = [ get_attribute_node(s_network, s_node, :cap) * 
                    sum(get_attribute_edge(s_network, get_edge(s_network, s_node, s_neighbor), :cap) for s_neighbor in neighbors(s_network, s_node)) 
                    for s_node in cluster]
            best_node = cluster[findmin(capacities_score)[2]]
            for v_node in v_subgraph.nodes_of_main_graph
                temporary_placement[v_node] = best_node
            end
        end
    

        for v_subgraph in vn_subgraphs

            s_subgraph = assignment_virtual_substrate_subgraphs[v_subgraph]
            sub_instance = Instance(v_subgraph.graph, s_subgraph.graph)

            
            # Additional cost thing
            additional_costs = []
            for v_node in vertices(v_subgraph.graph)
                current_addition_costs = zeros(nv(s_subgraph.graph))
                original_v_node = v_subgraph.nodes_of_main_graph[v_node]
                
                for v_edge in vn_decompo.v_edges_master
                    if src(v_edge) == original_v_node
                        placement_of_dst_node = temporary_placement[dst(v_edge)]
                        for s_node in vertices(s_subgraph.graph)
                            original_s_node = s_subgraph.nodes_of_main_graph[s_node]
                            current_addition_costs[s_node] += base_shortest_paths.dists[original_s_node, placement_of_dst_node] / 2
                        end
                    end
                    if dst(v_edge) == original_v_node
                        placement_of_dst_node = temporary_placement[src(v_edge)]
                        for s_node in vertices(s_subgraph.graph)
                            original_s_node = s_subgraph.nodes_of_main_graph[s_node]
                            current_addition_costs[s_node] += base_shortest_paths.dists[original_s_node, placement_of_dst_node] / 2
                        end
                    end
                end
                
    
                push!(additional_costs, current_addition_costs)
            end
    
            # GETTING THE SUBMAPPING
            if solver == "mepso"
                sub_mapping, cost = solve_mepso_custom(sub_instance, additional_costs; nb_particle=25, nb_iter=25, time_max=1., print_things=false)
            elseif solver == "local-search"
                result = solve_local_search_routing(sub_instance, additional_costs; nb_particle=1, nb_local_search=300)
                sub_mapping = result[:mapping]
                cost = result[:mapping_cost]
            elseif solver == "milp"
                result = solve_compact_addition_routing_cost(sub_instance, additional_costs; time_solver = 60)
                sub_mapping = result["mapping"]
                cost = result["solution_value"]
            else
                println("I don't know your solver. Using mepso.")
                sub_mapping, cost = solve_mepso_custom(sub_instance, additional_costs; nb_particle=25, nb_iter=25, time_max=1., print_things=false)
            end

            if isnothing(sub_mapping) # invalid submapping!
                print("A submapping failed. ")
                continue
            end
        
            
            real_cost = 0

            node_placement = []
            for v_node in vertices(v_subgraph.graph)
                original_s_node = s_subgraph.nodes_of_main_graph[sub_mapping.node_placement][v_node]
                append!(node_placement, original_s_node)

                original_v_node = v_subgraph.nodes_of_main_graph[v_node]
                temporary_placement[original_v_node] = original_s_node
    
                real_cost += s_network[original_s_node][:cost]

            end


            edge_routing = Dict()
            for v_edge in edges(v_subgraph.graph)
                used_edges = []
                for s_edge in sub_mapping.edge_routing[v_edge].edges
                    real_s_edge = get_edge(s_network_dir, s_subgraph.nodes_of_main_graph[src(s_edge)], s_subgraph.nodes_of_main_graph[dst(s_edge)])
                    push!(used_edges, real_s_edge)
                    real_cost += s_network_dir[src(real_s_edge), dst(real_s_edge)][:cost]
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


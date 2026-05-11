# This pricer aim at giving many columns, on a lot of places of the substrate network, quickly.
# it's kinda the "pave the graph" problematic.

# Todo : get several columns, and no need to get the optimal !
# Todo : update the pricer with dual costs !

includet("../../../utils/import_utils.jl")
includet("../../../utils/kahip_wrapper.jl")

using CPLEX, JuMP



struct PricerSubSubstrate
    nodes_of_original_graph
    vn_subgraph
    subinstance
    original_instance
    model
end

function Base.show(io::IO, pricer::PricerSubSubstrate)
    println(io, "Some pricer for subvn $(pricer.vn_subgraph.graph[][:name]) on subsn $(pricer.subinstance.s_network[][:name])")
end



function set_up_pricer_sn_decompo(instance, vn_subgraph, clusters)

    s_network = instance.s_network
    pricers = []
    for (i_cluster, cluster) in enumerate(clusters)
        sub_s_network = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")
        subinstance = Instance(vn_subgraph.graph, sub_s_network)

        model = Model(CPLEX.Optimizer)
        set_up_subpb(model, subinstance, instance, vn_subgraph, cluster)

        push!(pricers,  PricerSubSubstrate(cluster, vn_subgraph, subinstance, instance, model))
    end

    return pricers
end



# computing the sn decompo. It depends on the subvirtual network...
function get_sn_decompo(s_network, nb_clusters, nb_nodes_per_clusters)


    # 1 : Partitionner
    inbalance = 0.10
    partition = partition_kahip(s_network.graph, nb_clusters, inbalance)
    clusters = [Vector{Int64}() for i in 1:nb_clusters]
    for s_node in vertices(s_network)
        push!(clusters[partition[s_node]], s_node)
    end

    # 2 : Corriger

    # 2a) correction by removing unconnected sets
    for cluster in clusters
        simple_subgraph, vmap = induced_subgraph(s_network.graph, cluster)
        if !is_connected(simple_subgraph)
            #println("It's not connected ! :(")
            #println("At the beginning the cluster is : $cluster")
            components = connected_components(simple_subgraph)
            component_sorted = sort(components, by=x->length(x), rev=true)
            new_cluster = [vmap[i] for i in component_sorted[1]]
            for subcluster in component_sorted[2:length(component_sorted)]
                #Let's add all those nodes to a (most) connected subgraph
                nodes_original = [vmap[node] for node in subcluster]
                subgraph_neighbors = zeros(Int, nb_clusters)
                for node in nodes_original
                    for neighbor in neighbors(s_network, node)
                        if neighbor ∉ cluster
                            subgraph_neighbors[partition[neighbor]] += 1
                        end
                    end
                end
                most_connected_subgraph = sortperm(subgraph_neighbors, rev=true)
                append!(clusters[most_connected_subgraph[1]], nodes_original)
                #println("Well let's add $nodes_original to cluster $(clusters[most_connected_subgraph[1]])")
                filter!(e->e∉nodes_original, cluster)
            end
        end
    end

    # 2b) adding nodes. For now, any adjacent nodes will do.

    i_cluster = 1
    subgraphs = []


    for cluster in clusters

        all_neighbors = Dict()
        for s_node in cluster
            for neigh in neighbors(s_network, s_node)
                if neigh ∉ cluster
                    if neigh ∉ keys(all_neighbors)
                        all_neighbors[neigh] = 1
                    else
                        all_neighbors[neigh] += 1
                    end
                end
            end
        end

        added = []

        nb_nodes_with_capacity=0
        for s_node in cluster
            if s_network[s_node][:cap] >= 1
                nb_nodes_with_capacity+=1
            end
        end

        while nb_nodes_with_capacity < nb_nodes_per_clusters

            # ranking the neighbors
            ranking = sort(collect(keys(all_neighbors)), by = x->all_neighbors[x], rev=true)
            # add the most connected neighbor
            new_s_node = ranking[1]
            if s_network[new_s_node][:cap] >= 1
                nb_nodes_with_capacity+=1
            end
            push!(cluster, new_s_node)
            push!(added, new_s_node)
            delete!(all_neighbors, new_s_node)

            for neigh in neighbors(s_network, new_s_node)
                if neigh ∉ cluster
                    if neigh ∉ keys(all_neighbors)
                        all_neighbors[neigh] = 1
                    else
                        all_neighbors[neigh] += 1
                    end
                end
            end
        end

        #println("Well I have $(length(cluster)) whereas I need $(nb_nodes_per_clusters)")
        sub_s_network = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")
        i_cluster += 1
        push!(subgraphs, sub_s_network)

        #write_added_nodes(s_network.graph, cluster, added, sub_s_network[][:name])
    end



    return clusters

end



function set_up_subpb(model, subinstance, original_instance, v_subgraph, s_network_original_nodes)

    set_silent(model)
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Tolerances_MIPGap", 0.15)

    v_network = subinstance.v_network
    s_network_dir = subinstance.s_network_dir
    s_network = subinstance.s_network

    original_v_network = original_instance.v_network
    original_s_network = original_instance.s_network

    ### Variables
    @variable(model, x[vertices(v_network), vertices(s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);

    

    ### Objective
    placement_cost = @expression(model, sum( s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);




    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
    end

    # one to one : one virtual node per substrate node
    for s_node in vertices(s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end

    # Node capacities
    for s_node in vertices(s_network)
        @constraint(model, sum(x[v_node, s_node] * v_network[v_node][:dem] for v_node in vertices(v_network)) <= s_network[s_node][:cap])
    end

    ## Edges 
    
    # edge capacity (undirected version !)
    for s_edge in edges(s_network)
        @constraint(model, 
            sum( v_network[src(v_edge), dst(v_edge)][:dem] * (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
                for v_edge in edges(v_network)) 
            <= 
            s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(s_network)
        for v_edge in edges(v_network)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end

    
    ## Departure constraints
    
    for s_node in vertices(s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end
    
    
    # Outgoing edges cap: pretty stupid but useful. Still valid...
    j = 0
    for v_node in vertices(v_network)
        for s_node in vertices(s_network)
            necessary_bw = degree(v_network, v_node)
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
                j+=1
            end 
        end
    end
    
    #outgoing: keep in mind the original v_network, and s_network, connectivity, to add some new constraints?
    i=0
    for v_node in vertices(v_network)
        necessary_bw = degree(original_v_network, v_subgraph.nodes_of_main_graph[v_node])
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(original_s_network, s_network_original_nodes[s_node], neighbor) for neighbor in neighbors(original_s_network, s_network_original_nodes[s_node])]
            available_bw = sum(original_s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
                i +=1
            end 
        end
    end

    #=Some weird ass constraint?
    for s_node in vertices(s_network)
        s_edges_incident = [get_edge(original_s_network, s_network_original_nodes[s_node], neighbor) for neighbor in neighbors(original_s_network, s_network_original_nodes[s_node])]
        available_bw = sum(original_s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)

        @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network, s_node) for v_edge in edges(v_network))
                                + sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network, s_node) for v_edge in edges(v_network)) 
                                + sum(x[v_node, s_node] * (degree(original_v_network, v_subgraph.nodes_of_main_graph[v_node]) - degree(v_network, v_node)) for v_node in vertices(v_network)) 
                                    <= available_bw)
    end
    =#
    
end



function update_pricer_sn_decompo(vn_decompo, pricer, dual_costs)

    model = pricer.model

    s_node_original_graph = pricer.nodes_of_original_graph
    sub_s_network_dir = pricer.subinstance.s_network_dir
    original_s_network = pricer.original_instance.s_network
    sub_s_network = pricer.subinstance.s_network
    v_subgraph = pricer.vn_subgraph

    ### Objective
    placement_cost = @expression(model, 
        sum( ( sub_s_network_dir[s_node][:cost] - dual_costs.capacity_s_node[s_node_original_graph[s_node]] ) * v_subgraph.graph[v_node][:dem] * model[:x][v_node, s_node] 
            for v_node in vertices(v_subgraph.graph) for s_node in vertices(sub_s_network_dir) ))

    routing_cost = @expression(model, sum( 
        ( sub_s_network[src(s_edge), dst(s_edge)][:cost] - dual_costs.capacity_s_edge[get_edge(original_s_network, s_node_original_graph[src(s_edge)], s_node_original_graph[dst(s_edge)])] ) 
        * v_subgraph.graph[src(v_edge), dst(v_edge)][:dem] * (model[:y][v_edge, get_edge(sub_s_network_dir, src(s_edge), dst(s_edge))] + model[:y][v_edge, get_edge(sub_s_network_dir, dst(s_edge), src(s_edge))])
                for v_edge in edges(v_subgraph.graph) for s_edge in edges(sub_s_network) ))


            
    # flow conservation
    flow_conservation_cost = AffExpr(0.)

    for s_node in vertices(sub_s_network)
        original_node = s_node_original_graph[s_node]
        for connecting_edge in vn_decompo.v_edges_master
            if v_subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[src(connecting_edge)][v_subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    -dual_costs.flow_conservation[connecting_edge][original_node] , 
                    model[:x][v_node_subgraph, s_node])
            end
            if v_subgraph ∈ keys(vn_decompo.v_nodes_assignment[dst(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[dst(connecting_edge)][v_subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    +dual_costs.flow_conservation[connecting_edge][original_node], 
                    model[:x][v_node_subgraph, s_node])
            end
        end
    end


    # departure !
    departure_costs = AffExpr(0.)
    for s_node in vertices(sub_s_network)
        for connecting_edge in vn_decompo.v_edges_master
            if v_subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[src(connecting_edge)][v_subgraph]
                add_to_expression!(
                    departure_costs, 
                    -dual_costs.departure[connecting_edge][s_node_original_graph[s_node]], 
                    model[:x][v_node_subgraph, s_node])
            end
        end
    end


    @objective(model, Min, 
            -dual_costs.convexity[v_subgraph]
            + placement_cost + routing_cost 
            + flow_conservation_cost 
            + departure_costs);


    return
end




function solve_pricers_sn_decompo(pricer; time_limit = 1000)

    original_s_network_dir = pricer.original_instance.s_network_dir
    s_nodes_original = pricer.nodes_of_original_graph
    vn_subgraph = pricer.vn_subgraph
    v_network = vn_subgraph.graph
    s_network_dir = pricer.subinstance.s_network_dir
    model = pricer.model

    set_time_limit_sec(model, time_limit)
    optimize!(model)
    status = termination_status(model)

    if status == MOI.FEASIBLE_POINT || status == MOI.OPTIMAL
        # Get the solution
        x_values = value.(model[:x])
        y_values = value.(model[:y])
        true_cost = 0.

        node_placement = []
        for v_node in vertices(v_network)
            for s_node in vertices(s_network_dir)
                if x_values[v_node, s_node] > 0.99
                    real_s_node = s_nodes_original[s_node]
                    append!(node_placement, real_s_node)
                    true_cost += original_s_network_dir[real_s_node][:cost]
                end
            end
        end


        edge_routing = Dict()
        for v_edge in edges(v_network)
            if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
                edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
            end
            used_edges = []
            for s_edge in edges(s_network_dir)
                if y_values[v_edge, s_edge] > 0.99
                    real_s_edge = get_edge(original_s_network_dir, s_nodes_original[src(s_edge)], s_nodes_original[dst(s_edge)])
                    push!(used_edges, real_s_edge)
                    true_cost += original_s_network_dir[src(real_s_edge), dst(real_s_edge)][:cost]
                end
            end
            edge_routing[v_edge] = order_path(original_s_network_dir, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
        end
        mapping = Mapping(v_network, original_s_network_dir, node_placement, edge_routing)
        reduced_cost = objective_value(model)
        return mapping, true_cost, reduced_cost
    else
        return nothing, 9999999, 999999
    end


end




function solve_pricers_sn_decompo_several_sols(pricer; time_limit = 100)

    original_s_network_dir = pricer.original_instance.s_network_dir
    s_nodes_original = pricer.nodes_of_original_graph
    vn_subgraph = pricer.vn_subgraph
    v_network = vn_subgraph.graph
    s_network_dir = pricer.subinstance.s_network_dir
    model = pricer.model

    set_time_limit_sec(model, time_limit)
    optimize!(model)
    status = termination_status(model)

    if status == MOI.FEASIBLE_POINT || status == MOI.OPTIMAL
        # Get the solution
        x_values = value.(model[:x])
        y_values = value.(model[:y])
        true_cost = 0.

        node_placement = []
        for v_node in vertices(v_network)
            for s_node in vertices(s_network_dir)
                if x_values[v_node, s_node] > 0.99
                    real_s_node = s_nodes_original[s_node]
                    append!(node_placement, real_s_node)
                    true_cost += original_s_network_dir[real_s_node][:cost]
                end
            end
        end


        edge_routing = Dict()
        for v_edge in edges(v_network)
            if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
                edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
            end
            used_edges = []
            for s_edge in edges(s_network_dir)
                if y_values[v_edge, s_edge] > 0.99
                    real_s_edge = get_edge(original_s_network_dir, s_nodes_original[src(s_edge)], s_nodes_original[dst(s_edge)])
                    push!(used_edges, real_s_edge)
                    true_cost += original_s_network_dir[src(real_s_edge), dst(real_s_edge)][:cost]
                end
            end
            edge_routing[v_edge] = order_path(original_s_network_dir, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
        end
        mapping = Mapping(v_network, original_s_network_dir, node_placement, edge_routing)
        reduced_cost = objective_value(model)
        return mapping, true_cost, reduced_cost
    else
        return nothing, 9999999, 999999
    end


end




includet("../../../utils/import_utils.jl")




# We get the central nodes of the network, and then we look how all the nodes are far from these nodes. 
# Depending on this, we compute a placement that should be decent, as it should be compact 

function get_decent_placements(instance)
    central_s_node = argmax(betweenness_centrality(instance.s_network))
    central_v_nodes = Dict() 
    for v_network in instance.v_networks
        central_v_nodes[v_network] = argmax(betweenness_centrality(v_network))
    end

    shortest_path_substrate = get_shortest_paths(instance.s_network, 1);
    shortest_path_virtuals = Dict()
    for v_network in instance.v_networks
        shortest_path_virtuals[v_network] = get_shortest_paths(v_network, 1)
    end

    penalty_substrate = []
    for s_node in vertices(instance.s_network)
        shortest_path = shortest_path_substrate[(central_s_node, s_node)]
        if length(shortest_path)>0
            push!(penalty_substrate, shortest_path[1].cost)
        else
            push!(penalty_substrate, 9999) #it should not happen, but so far it does
        end
    end

    distance_virtuals = Dict()
    for v_network in instance.v_networks
        distance_virtual = []
        for v_node in vertices(v_network)
            shortest_path = shortest_path_substrate[(central_v_nodes[v_network], v_node)]
            if length(shortest_path)>0
                if shortest_path[1].cost == 0
                    push!(distance_virtual, 1)
                else
                    push!(distance_virtual, shortest_path[1].cost)
                end
            else
                push!(distance_virtual, 1) #it should not happen, but so far it does
            end
        end
        distance_virtuals[v_network] = distance_virtual
    end

    
    
    
    #### Model
    model_placement = Model(CPLEX.Optimizer)
    set_attribute(model_placement, "CPX_PARAM_EPINT", 1e-8)

    ### Variables
    x_variables = @variable(model_placement, x[v_network in instance.v_networks, vertices(v_network), vertices(instance.s_network)], binary=true);

    ### Objective
    placement_cost = @expression(model_placement, sum( 
        instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_network, v_node, s_node] 
        for v_network in instance.v_networks for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    @objective(model_placement, Min, placement_cost);

    ### Constraints

    # one substrate node per virtual node
    for v_network in instance.v_networks
        for v_node in vertices(v_network)
            @constraint(model_placement, sum(x[v_network, v_node, s_node] for s_node in vertices(instance.s_network)) == 1)
        end
    end

    # node capacity
    for s_node in vertices(instance.s_network)
        @constraint(model_placement, 
            sum( v_network[v_node][:dem] * x[v_network, v_node, s_node] 
                for v_network in instance.v_networks for v_node in vertices(v_network) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end

    # Solving
    set_time_limit_sec(model_placement, 30.)
    set_silent(model_placement)
    optimize!(model_placement)

    # Get the solution
    x_values = value.(x_variables);
    node_placements = Dict()
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            current_placement = []
            for s_node in vertices(instance.s_network)
                if x_values[v_network, v_node, s_node] > 0.99
                    push!(current_placement, 1)
                else
                    push!(current_placement, 0)
                end
            end
            push!(node_placement, current_placement)
        end
        node_placements[v_network] = node_placement
    end



    return node_placements

end
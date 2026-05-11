struct NetworkDecomposition
    node_partitionning
    subgraphs
    node_assignment
    connecting_edges
    mappings
end

struct MasterProblem
    model
    vn_decompos
    lambdas
end


struct DualCosts
    convexity
    capacity_s_node
    capacity_s_edge
    flow_conservation
    departure
end


struct MappingDecompoFrac
    v_network
    s_network
    mapping_selection
    mapping_costs
    connecting_edge_routing
    connecting_cost
end


function MappingDecompoFrac(v_network, s_network, mapping_selection, connecting_edge_routing)
    # TODO!!! :( 
    placement_costs = 0
    routing_costs = 0
    return MappingDecompoFrac(v_network, s_network, mapping_selection, placement_costs, connecting_edge_routing, routing_costs)
end

function Base.show(io::IO, mapping::MappingDecompoFrac)

    println(io, "Overall mapping of subgraphs costs : " * string(mapping.mapping_costs))
    for i_subgraph in 1:length(mapping.mapping_selection)
        println(io, "Mapping of subgraph " * string(i_subgraph))
        i_map = 1
        for mapping_subgraph in keys(mapping.mapping_selection[i_subgraph])
            println("Mapping nb " * string(i_map) * " with a value " * string( round(mapping.mapping_selection[i_subgraph][mapping_subgraph] * 1000) / 1000) )
            println(mapping_subgraph)
            println(mapping.mapping_selection[i_subgraph][mapping_subgraph])
        end
    end

    println("\n\n Edge routing : ")
    for v_edge in keys(mapping.connecting_edge_routing)
        println("For connecting edge " * string(v_edge) * " : ")
        for s_edge in keys(mapping.connecting_edge_routing[v_edge])
            println("The edge " * string(s_edge) * " : " * string(mapping.connecting_edge_routing[v_edge][s_edge]))
        end
    end
end

# utils
function my_induced_subgraph(meta_graph::MetaGraph, selector, name)
    induced_graph, osef = induced_subgraph(meta_graph.graph, selector)

    induced_metagraph = MetaGraph(
        DiGraph(),
        Int,
        Dict,
        Dict,
        Dict(:name => name, :type => meta_graph[][:type], :directed => true)
    )

    for node in vertices(induced_graph)
        add_vertex!(induced_metagraph, node, meta_graph[selector[node]])
    end

    for edge in edges(induced_graph)
        add_edge!(induced_metagraph, src(edge), dst(edge), meta_graph[selector[src(edge)], selector[dst(edge)]])
    end

    return induced_metagraph
end



function solve_integer_with_forbidden_nodes(instance, forbidden_nodes)

    #### Model
    model = Model(CPLEX.Optimizer)
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)

    ### Variables
    @variable(model, x[v_network in instance.v_networks, vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[v_network in instance.v_networks, edges(v_network), edges(instance.s_network)], binary=true);

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_network, v_node, s_node] 
        for v_network in instance.v_networks for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( instance.s_network[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge]
        for v_network in instance.v_networks for v_edge in edges(v_network) for s_edge in edges(instance.s_network) ))
    @objective(model, Min, placement_cost + routing_cost);


    ### Constraints     

    one_to_one = true
    departure_cst = true

    ## Nodes

    # one substrate node per virtual node
    for v_network in instance.v_networks
        for v_node in vertices(v_network)
            @constraint(model, sum(x[v_network, v_node, s_node] for s_node in vertices(instance.s_network)) == 1)
        end
    end

    # if one to one : one virtual node per substrate node
    if one_to_one
        for s_node in vertices(instance.s_network)
            for v_network in instance.v_networks
                @constraint(model, sum(x[v_network, v_node, s_node] for v_node in vertices(v_network)) <= 1)
            end
        end
    end

    # node capacity
    for s_node in vertices(instance.s_network)
        @constraint(model, 
            sum( v_network[v_node][:dem] * x[v_network, v_node, s_node] 
                for v_network in instance.v_networks for v_node in vertices(v_network) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end


    ## Edges 
    
    # edge capacity
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge] 
                for v_network in instance.v_networks for v_edge in edges(v_network)) 
            <= 
            instance.s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(instance.s_network)
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                @constraint(model, 
                    x[v_network, src(v_edge), s_node] - x[v_network, dst(v_edge), s_node] 
                    <=
                    sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) - 
                        sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node))
                )
            end
        end
    end


    ## Additional constraints : Node + Edge
    if one_to_one
        if departure_cst
            for s_node in vertices(instance.s_network)
                for v_network in instance.v_networks
                    for v_node in vertices(v_network)
                        for v_edge in get_out_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                        end
                        for v_edge in get_in_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                        end
                    end
                end
            end
        end
    end

    ## Additional constraints : Undirected edges, thus the value both way is the same
    undirected = true
    if undirected
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                for s_edge in edges(instance.s_network)
                    @constraint(model, y[v_network, v_edge, s_edge] == y[v_network, get_edge(v_network, dst(v_edge), src(v_edge)), 
                                                get_edge(instance.s_network, dst(s_edge), src(s_edge))])
                end
            end
        end
    end

    # Additional constraint: the forbidden nodes can not host any virtual nodes
    for s_node in forbidden_nodes
        for v_network in instance.v_networks
            for v_node in vertices(v_network)
                @constraint(model, x[v_network, v_node, s_node] == 0)
            end
        end
    end

    # Solving
    set_silent(model)
    optimize!(model)

    # Get the solution
    x_values = value.(model[:x])
    y_values = value.(model[:y])
    mappings = []
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            for s_node in vertices(instance.s_network)
                if x_values[v_network, v_node, s_node] > 0.99
                    append!(node_placement, s_node)
                end
            end
        end

        edge_routing = Dict()
        for v_edge in edges(v_network)
            if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
                edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
            end
            used_edges = []
            for s_edge in edges(instance.s_network)
                if y_values[v_network, v_edge, s_edge] > 0.99
                    push!(used_edges, s_edge)
                end
            end
            edge_routing[v_edge] = order_path(instance.s_network, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
        end
        m = Mapping(v_network, instance.s_network, node_placement, edge_routing)
        push!(mappings, m)
    end

    return mappings
end


# for unit capacity only (need small adaptation, remove the used cap)
function get_initial_set_of_mappings(instance, vn_decompos)
    forbidden_nodes = []
    
    for v_network in instance.v_networks

        for subgraph in vn_decompos[v_network].subgraphs
            current_instance = InstanceVNE([subgraph], instance.s_network)
            mapping = solve_integer_with_forbidden_nodes(current_instance, forbidden_nodes)[1]
            for s_node in mapping.node_placement
                push!(forbidden_nodes, s_node)
            end
            push!(vn_decompos[v_network].mappings[subgraph], mapping);
        end
    end
end



### A better initialization: trying to put each vn a bit everywhere
# does not work for overlapping decompo for now
function get_initial_set_of_mappings_better(instance, vn_decompos)

    println("YO LET'S COOK MORE COLUMNS")

    nb_column_per_subgraph = 10

    nb_s_nodes = length(vertices(instance.s_network))

    #substrate_nodes = shuffle(L) # for random
    substrate_nodes = 1:nb_s_nodes

    # Compute the approximate size of each partition
    base_size = div(nb_s_nodes, nb_column_per_subgraph)
    extra = nb_s_nodes % nb_column_per_subgraph

    # Partition the shuffled list into roughly equal parts
    groups = []
    start_idx = 1
    for i in 1:nb_column_per_subgraph
        # Distribute the remainder across the first few groups
        group_size = base_size + (i <= extra ? 1 : 0)
        end_idx = start_idx + group_size - 1
        push!(groups, substrate_nodes[start_idx:end_idx])
        start_idx = end_idx + 1
    end

    #println("groups: $groups")
    # get the most central node
    # assign it somewhere on the selected nodes.
    # => nodes at random ? It would be better through some gentle clustering ?
    # No need to be just one node, it should be much better if it has some choice. I think random would be good to begin with. 
    # solve the plne, add the column.



    for v_network in instance.v_networks
        for subgraph in vn_decompos[v_network].subgraphs
            current_instance = InstanceVNE([subgraph], instance.s_network)
            println("Let's have fun !")
            for group in groups
                # the virtual node should be not too connected to make it easy ?
                print(group)
                placement_restriction = Dict()
                placement_restriction[1] = group
                # create the plne
               

                mapping = solve_integer_with_placement_restriction(current_instance, placement_restriction)
                
                if mapping !== nothing
                    push!(vn_decompos[v_network].mappings[subgraph], mapping);
                end
            end
        end
    end

end

# Todo: use the base function...
function solve_integer_with_placement_restriction(instance, placement_restriction)

    print("En zarbi")
    #### Model
    model = Model(CPLEX.Optimizer)

    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)

    ### Variables
    @variable(model, x[v_network in instance.v_networks, vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[v_network in instance.v_networks, edges(v_network), edges(instance.s_network)], binary=true);


    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_network, v_node, s_node] 
        for v_network in instance.v_networks for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))

    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_network, v_node, s_node]
        for v_network in instance.v_networks for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( instance.s_network[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge]
        for v_network in instance.v_networks for v_edge in edges(v_network) for s_edge in edges(instance.s_network) ))
    @objective(model, Min, placement_cost + routing_cost);

    println("How we lookin")

    ### Constraints     

    one_to_one = true
    departure_cst = true


    ### ======== Additional constraints
    ### placement_restriction
    for v_network in instance.v_networks
        for v_node in keys(placement_restriction)
            @constraint(model, sum(x[v_network, v_node, s_node] for s_node in placement_restriction[v_node]) == 1)
        end
    end


    ### =========== Nodes constraints

    # one substrate node per virtual node
    for v_network in instance.v_networks
        for v_node in vertices(v_network)
            @constraint(model, sum(x[v_network, v_node, s_node] for s_node in vertices(instance.s_network)) == 1)
        end
    end

    # if one to one : one virtual node per substrate node
    if one_to_one
        for s_node in vertices(instance.s_network)
            for v_network in instance.v_networks
                @constraint(model, sum(x[v_network, v_node, s_node] for v_node in vertices(v_network)) <= 1)
            end
        end
    end

    # node capacity
    for s_node in vertices(instance.s_network)
        @constraint(model, 
            sum( v_network[v_node][:dem] * x[v_network, v_node, s_node] 
                for v_network in instance.v_networks for v_node in vertices(v_network) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end


    ### ========== Edges constraints 
    
    # edge capacity
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge] 
                for v_network in instance.v_networks for v_edge in edges(v_network)) 
            <= 
            instance.s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(instance.s_network)
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                @constraint(model, 
                    x[v_network, src(v_edge), s_node] - x[v_network, dst(v_edge), s_node] 
                    <=
                    sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) - 
                        sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node))
                )
            end
        end
    end


    ### ========= Additional constraints
    
    ## Departure constraint
    if one_to_one
        if departure_cst
            for s_node in vertices(instance.s_network)
                for v_network in instance.v_networks
                    for v_node in vertices(v_network)
                        for v_edge in get_out_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                        end
                        for v_edge in get_in_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                        end
                    end
                end
            end
        end
    end

    ## Symmetric edges routings (because undirected !)
    undirected = true
    if undirected
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                for s_edge in edges(instance.s_network)
                    @constraint(model, y[v_network, v_edge, s_edge] == y[v_network, get_edge(v_network, dst(v_edge), src(v_edge)), 
                                                get_edge(instance.s_network, dst(s_edge), src(s_edge))])
                end
            end
        end
    end

    println("How we lookin")
    # Solving
    set_silent(model)
    optimize!(model)

    status = termination_status(model)

    if status != MOI.OPTIMAL
        println("Infeasible or unfinished: $status")
        return
    end

    print("Bah alors ?")
    # Get the solution
    x_values = value.(model[:x])
    y_values = value.(model[:y])
    columns = []
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            for s_node in vertices(instance.s_network)
                if x_values[v_network, v_node, s_node] > 0.99
                    append!(node_placement, s_node)
                end
            end
        end

        edge_routing = Dict()
        for v_edge in edges(v_network)
            if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
                edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
            end
            used_edges = []
            for s_edge in edges(instance.s_network)
                if y_values[v_network, v_edge, s_edge] > 0.99
                    push!(used_edges, s_edge)
                end
            end
            edge_routing[v_edge] = order_path(instance.s_network, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
        end
        m = Mapping(v_network, instance.s_network, node_placement, edge_routing)
        println(m)

        return m

    end

end





############ MASTER PROBLEM

function set_up_master_problem(instance, vn_decompos)
    model = Model(CPLEX.Optimizer)
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)
    ### Variables
    
    
    @variable(model, y[
        v_network in instance.v_networks, 
        v_edge in vn_decompos[v_network].connecting_edges, 
        s_edge in edges(instance.s_network)], 
        binary=true);
    
    ### Objective
    connecting_cost = @expression(model, sum( instance.s_network[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge]
        for v_network in instance.v_networks for v_edge in vn_decompos[v_network].connecting_edges for s_edge in edges(instance.s_network) ))
    
    @objective(model, Min, connecting_cost);

    ### Constraints

    # convexity constraints
    @constraint(
        model, 
        mapping_selec[v_network in instance.v_networks, subgraph in vn_decompos[v_network].subgraphs],
        0 == 1
    );


    # capacity of substrate noeuds !
    @constraint(
        model,
        capacity_s_node[s_node in vertices(instance.s_network)],
        0 <= instance.s_network[s_node][:cap]
    );


    # capacity of substrate edges
    @constraint(
        model,
        capacity_s_edge[s_edge in edges(instance.s_network)],
        sum( v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge] 
            for v_network in instance.v_networks for v_edge in vn_decompos[v_network].connecting_edges)
        + 0
        <= instance.s_network[src(s_edge), dst(s_edge)][:cap]
    );


    # flow conservation constraints
    @constraint(
        model,
        flow_conservation[v_network in instance.v_networks, connecting_edge in vn_decompos[v_network].connecting_edges, s_node in vertices(instance.s_network)],
        0 == 
        sum( y[v_network, connecting_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node))
        - sum( y[v_network, connecting_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node))
    );


    # Departure constraints
    @constraint(
        model, 
        departure[v_network in instance.v_networks, v_edge in vn_decompos[v_network].connecting_edges, s_node in vertices(instance.s_network)],
        sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) >=
        0
    )

    ## Additional constraints : Undirected edges, thus the value both way is the same
    undirected = true
    if undirected
        for v_network in instance.v_networks
            for v_edge in vn_decompos[v_network].connecting_edges
                for s_edge in edges(instance.s_network)
                    @constraint(model, y[v_network, v_edge, s_edge] == y[v_network, get_edge(v_network, dst(v_edge), src(v_edge)), 
                                                get_edge(instance.s_network, dst(s_edge), src(s_edge))])
                end
            end
        end
    end

    lambdas = Dict()
    for v_network in instance.v_networks
        lambdas[v_network] = Dict()
        for subgraph in vn_decompos[v_network].subgraphs
            lambdas[v_network][subgraph] = Dict()
        end
    end
    return MasterProblem(model, vn_decompos, lambdas)
end


function add_column(master_problem, instance, v_network, subgraph, mapping)

    lambda = @variable(master_problem.model, binary=true);

    master_problem.lambdas[v_network][subgraph][mapping] = lambda

    set_objective_coefficient(master_problem.model, lambda, mapping.node_placement_cost + mapping.edge_routing_cost)

    # convexity
    set_normalized_coefficient(master_problem.model[:mapping_selec][v_network, subgraph], lambda, 1)

    
    # capacities
    for s_node in vertices(instance.s_network)
        set_normalized_coefficient(master_problem.model[:capacity_s_node][s_node], lambda, mapping.s_node_usage[s_node])
    end
    for s_edge in edges(instance.s_network)
        set_normalized_coefficient(master_problem.model[:capacity_s_edge][s_edge], lambda, mapping.s_edge_usage[s_edge])
    end

    
    # flow conservation
    for connecting_edge in master_problem.vn_decompos[v_network].connecting_edges
        if subgraph == master_problem.vn_decompos[v_network].node_assignment[src(connecting_edge)][1]
            set_normalized_coefficient(master_problem.model[:flow_conservation][v_network, connecting_edge, mapping.node_placement[master_problem.vn_decompos[v_network].node_assignment[src(connecting_edge)][2]]], 
                lambda, 
                1)
        end
        if subgraph == master_problem.vn_decompos[v_network].node_assignment[dst(connecting_edge)][1]
            set_normalized_coefficient(master_problem.model[:flow_conservation][v_network, connecting_edge, mapping.node_placement[master_problem.vn_decompos[v_network].node_assignment[dst(connecting_edge)][2]]], 
                lambda,
                -1)
        end
    end
    
    # departure
    for connecting_edge in master_problem.vn_decompos[v_network].connecting_edges
        if subgraph == master_problem.vn_decompos[v_network].node_assignment[src(connecting_edge)][1]
            set_normalized_coefficient(master_problem.model[:departure][v_network, connecting_edge, mapping.node_placement[master_problem.vn_decompos[v_network].node_assignment[src(connecting_edge)][2]]], 
                lambda, 
                -1)
        end
    end


end


######## SUBPROBLEMS




function get_duals(instance, vn_decompos, master_problem)
    
    convexity = Dict()
    flow_conservation = Dict()
    departure = Dict()

    for v_network in instance.v_networks

        convexity[v_network] = Dict()
        for subgraph in vn_decompos[v_network].subgraphs
            convexity[v_network][subgraph] = dual(master_problem.model[:mapping_selec][v_network, subgraph])
        end

        flow_conservation[v_network] = Dict()
        for v_edge in vn_decompos[v_network].connecting_edges
            flow_conservation[v_network][v_edge] = Dict()
            for s_node in vertices(instance.s_network)
                flow_conservation[v_network][v_edge][s_node] = dual(master_problem.model[:flow_conservation][v_network, v_edge, s_node])
            end
        end

        departure[v_network] = Dict()
        for v_edge in vn_decompos[v_network].connecting_edges
            departure[v_network][v_edge] = Dict()
            for s_node in vertices(instance.s_network)
                departure[v_network][v_edge][s_node] = dual(master_problem.model[:departure][v_network, v_edge, s_node])
            end
        end


    end
    
    capacity_s_node = Dict()

    for s_node in vertices(instance.s_network)
        capacity_s_node[s_node]  = dual(master_problem.model[:capacity_s_node][s_node])
    end

    capacity_s_edge = Dict()
    for s_edge in edges(instance.s_network)
        capacity_s_edge[s_edge]  = dual(master_problem.model[:capacity_s_edge][s_edge])
    end


    return DualCosts(convexity, capacity_s_node, capacity_s_edge, flow_conservation, departure)
end


struct SubProblem
    model
    vn_decompos
    s_network
    v_network
    subgraph
end


function set_up_pricer(instance, vn_decompos, v_network_used, subgraph)

    s_network = instance.s_network

    #### Model
    model = Model(CPLEX.Optimizer)
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)

    ### Variables
    @variable(model, x[v_node in vertices(subgraph), s_node in vertices(s_network)], binary=true);
    @variable(model, y[v_edge in edges(subgraph), s_edge in edges(s_network)], binary=true);


    ### Constraints


    one_to_one = true
    departure_cst = true

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(subgraph)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
    end

    # if one to one : one virtual node per substrate node
    if one_to_one
        for s_node in vertices(s_network)
            @constraint(model, sum(x[v_node, s_node] for v_node in vertices(subgraph)) <= 1)
        end
    end



    # node capacity
    for s_node in vertices(s_network)
        @constraint(model, 
            sum( subgraph[v_node][:dem] * x[v_node, s_node] 
                for v_node in vertices(subgraph) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end


    ## Edges 
    
    # edge capacity
    for s_edge in edges(s_network)
        @constraint(model, 
            sum( subgraph[src(v_edge), dst(v_edge)][:dem] * y[v_edge, s_edge] 
                for v_edge in edges(subgraph)) 
            <= 
            s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(s_network)
        for v_edge in edges(subgraph)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                <=
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network, s_node))
            )
        end
    end


    ## Additional constraints : Node + Edge
    if one_to_one && departure_cst
        for s_node in vertices(s_network)
            for v_node in vertices(subgraph)
                for v_edge in get_out_edges(subgraph, v_node)
                    @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network, s_node)) >= x[v_node, s_node])
                end
                for v_edge in get_in_edges(subgraph, v_node)
                    @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network, s_node)) >= x[v_node, s_node])
                end
            end
        end
    end

    ## Additional constraints : Undirected edges, thus the value both way is the same
    undirected = true
    if undirected
        for v_edge in edges(subgraph)
            for s_edge in edges(s_network)
                @constraint(model, y[v_edge, s_edge] == y[get_edge(subgraph, dst(v_edge), src(v_edge)), 
                                            get_edge(s_network, dst(s_edge), src(s_edge))])
            end
        end
    end



    return SubProblem(model, vn_decompos, instance.s_network, v_network_used, subgraph);
end


function update_solve_pricer(instance, pricer, dual_costs)

    model = pricer.model
    subgraph = pricer.subgraph
    vn_decompos = pricer.vn_decompos

    ### Objective
    placement_cost = @expression(model, sum( ( pricer.s_network[s_node][:cost] - dual_costs.capacity_s_node[s_node] ) * subgraph[v_node][:dem] * model[:x][v_node, s_node] 
        for v_node in vertices(subgraph) for s_node in vertices(pricer.s_network) ))

    routing_cost = @expression(model, sum( ( pricer.s_network[src(s_edge), dst(s_edge)][:cost] - dual_costs.capacity_s_edge[s_edge] ) * subgraph[src(v_edge), dst(v_edge)][:dem] * model[:y][v_edge, s_edge] 
        for v_edge in edges(subgraph) for s_edge in edges(pricer.s_network) ))


            
    # flow conservation
    flow_conservation_cost = AffExpr(0.)


    for s_node in vertices(pricer.s_network)
        for connecting_edge in vn_decompos[pricer.v_network].connecting_edges
            if subgraph == vn_decompos[pricer.v_network].node_assignment[src(connecting_edge)][1]
                add_to_expression!(flow_conservation_cost, -dual_costs.flow_conservation[pricer.v_network][connecting_edge][s_node], model[:x][vn_decompos[pricer.v_network].node_assignment[src(connecting_edge)][2], s_node])
            end
            if subgraph == vn_decompos[pricer.v_network].node_assignment[dst(connecting_edge)][1]
                add_to_expression!(flow_conservation_cost, +dual_costs.flow_conservation[pricer.v_network][connecting_edge][s_node], model[:x][vn_decompos[pricer.v_network].node_assignment[dst(connecting_edge)][2], s_node])
            end
        end
    end
            

    departure_costs = AffExpr(0.)
    for s_node in vertices(pricer.s_network)
        for connecting_edge in vn_decompos[pricer.v_network].connecting_edges
            if subgraph == vn_decompos[pricer.v_network].node_assignment[src(connecting_edge)][1]
                add_to_expression!(flow_conservation_cost, dual_costs.departure[pricer.v_network][connecting_edge][s_node], model[:x][vn_decompos[pricer.v_network].node_assignment[src(connecting_edge)][2], s_node])
            end
        end
    end


    @objective(model, Min, -dual_costs.convexity[pricer.v_network][subgraph] + placement_cost + routing_cost + flow_conservation_cost + departure_costs);


    set_silent(model)
    optimize!(model)


    # Get the solution
    x_values = value.(pricer.model[:x])
    y_values = value.(pricer.model[:y])
    node_placement = []
    for v_node in vertices(subgraph)
        for s_node in vertices(pricer.s_network)
            if x_values[v_node, s_node] > 0.99
                append!(node_placement, s_node)
            end
        end
    end

    edge_routing = Dict()
    for v_edge in edges(subgraph)
        if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
            edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
        end
        used_edges = []
        for s_edge in edges(instance.s_network)
            if y_values[v_edge, s_edge] > 0.99
                push!(used_edges, s_edge)
            end
        end
        edge_routing[v_edge] = order_path(pricer.s_network, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
    end
    mapping = Mapping(subgraph, instance.s_network, node_placement, edge_routing)
    
    dual_value = objective_value(model)
        
    return mapping, dual_value


end



function vn_decompo(instance, node_partitionning)

    vn_decompos = Dict()
    
    println("Starting solving...")

    for i_vn in 1:length(instance.v_networks)
        
        node_assignment = Dict()
        connecting_edges = []    
        subgraphs = []

        for i_subgraph in 1:length(node_partitionning[i_vn])
            nodes = node_partitionning[i_vn][i_subgraph]
            subgraph = my_induced_subgraph(instance.v_networks[i_vn], nodes, "subgraph_" * string(i_subgraph))
            for i_node in 1:length(nodes)
                node_assignment[nodes[i_node]] = [subgraph, i_node]
            end
            print(subgraph)
            push!(subgraphs, subgraph)
        end

        for v_edge in edges(instance.v_networks[1])
            if node_assignment[src(v_edge)][1] != node_assignment[dst(v_edge)][1]
                push!(connecting_edges, v_edge)
            end
        end
    
        mappings = Dict()
        for subgraph in subgraphs
            mappings[subgraph] = []
        end

        vn_decompos[instance.v_networks[i_vn]] = NetworkDecomposition(node_partitionning[i_vn], subgraphs, node_assignment, connecting_edges, mappings)

    end

    master_problem = set_up_master_problem(instance, vn_decompos)
    print("Master problem set... ")
    nb_columns = 0

    get_initial_set_of_mappings(instance, vn_decompos)
    get_initial_set_of_mappings_better(instance, vn_decompos)
    for v_network in instance.v_networks
        for subgraph in vn_decompos[v_network].subgraphs
            for mapping in vn_decompos[v_network].mappings[subgraph]
                add_column(master_problem, instance, v_network, subgraph, mapping)
                nb_columns += 1
                println("Column of price $(mapping.node_placement_cost + mapping.edge_routing_cost)")

            end
        end
    end
    print("Initial set of columns generated... ")

    pricers = Dict()
    for v_network in instance.v_networks
        pricers[v_network] = Dict()
        for subgraph in vn_decompos[v_network].subgraphs
            pricers[v_network][subgraph] = set_up_pricer(instance, vn_decompos, v_network, subgraph)
        end
    end


    # GENERATION DE COLONNES !
    set_silent(master_problem.model)
    time_master = 0
    time_subproblems = 0
    keep_on = true
    nb_iter = 0
    best_LG = -100000
    unrelax = nothing
    print("Starting column generation: \n")
    while keep_on && nb_iter < 100

        print("Iter " * string(nb_iter))
        if nb_iter > 0
            unrelax()
        end

        time1 = time()
        unrelax = relax_integrality(master_problem.model)
        optimize!(master_problem.model)
        dual_costs = get_duals(instance, vn_decompos, master_problem)
        time2 = time()

        nb_iter += 1
        keep_on = false

        print(", CG value : " * string( floor(objective_value(master_problem.model)* 1000) / 1000) )
        total_subpb_obj = 0
        for v_network in instance.v_networks
            for subgraph in vn_decompos[v_network].subgraphs
                mapping, obj_value = update_solve_pricer(instance, pricers[v_network][subgraph], dual_costs)
                total_subpb_obj += obj_value
                if obj_value < -0.0001
                    keep_on = true 
                    add_column(master_problem, instance, v_network, subgraph, mapping)
                    nb_columns += 1
                    push!(vn_decompos[v_network].mappings[subgraph], mapping)
                end
            end
        end
        time3 = time()
        

        # Calculating LG bound
        total_dual_value = 0
        for v_network in instance.v_networks
            for subgraph in vn_decompos[v_network].subgraphs
                total_dual_value += dual_costs.convexity[v_network][subgraph]
            end
        end
        for s_node in vertices(instance.s_network)
            total_dual_value += dual_costs.capacity_s_node[s_node] * instance.s_network[s_node][:cap]
        end
        for s_edge in edges(instance.s_network)
            total_dual_value += dual_costs.capacity_s_edge[s_edge] * instance.s_network[src(s_edge), dst(s_edge)][:cap]
        end
        LG_value =  total_dual_value + total_subpb_obj
        if LG_value > best_LG
            best_LG = LG_value
        end
        print(", Lagrangian Bound: " * string(floor(best_LG * 1000 ) / 1000 ))


        print(", Nb Columns: " * string(nb_columns))
        time_master += time2 - time1
        time_subproblems += time3 - time2
        print("\n")
    end



    ##### RECUPERATION DES SOLUTIONS
    if nb_iter > 0
        unrelax()
    end
    unrelax = relax_integrality(master_problem.model)
    optimize!(master_problem.model)

    println("________________ \nCG finished\nFinal value: " * string(objective_value(master_problem.model)))
    println("Time in MP: " * string(time_master) * ", time in SP: " * string(time_subproblems))
    y_values = value.(master_problem.model[:y])

    mappings = []
    for v_network in instance.v_networks
        mapping_selec = []
        for subgraph in vn_decompos[v_network].subgraphs
            subgraph_map = Dict()
            for mapping in vn_decompos[v_network].mappings[subgraph]
                if value.(master_problem.lambdas[v_network][subgraph][mapping]) > 0.01
                    subgraph_map[mapping] = value.(master_problem.lambdas[v_network][subgraph][mapping])
                end
            end
            push!(mapping_selec, subgraph_map)
        end
        
        edge_routing = Dict()
        for v_edge in vn_decompos[v_network].connecting_edges
            edges_used = Dict()
            for s_edge in edges(instance.s_network)
                if y_values[v_network, v_edge, s_edge] > 0.01
                    edges_used[s_edge] = y_values[v_network, v_edge, s_edge]
                end
            end
            edge_routing[v_edge] = edges_used
        end
        
        mapping = MappingDecompoFrac(v_network, instance.s_network, mapping_selec, edge_routing);
        push!(mappings, mapping)
    end


    ####### RESOLUTION ENTIERE
    unrelax()
    optimize!(master_problem.model)
    println("Value integer : " * string(objective_value(master_problem.model)))


    return mappings

end







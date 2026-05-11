

using Revise, JuMP, CPLEX, Gurobi
includet("../../utils/import_utils.jl")


struct Compact_Formulation
    model
    x_variables
    y_variables
end



function set_up_problem(instance, model)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    ### Variables
    @variable(model, x[vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);

    

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);




    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(instance.s_network)) == 1)
    end

    # one to one : one virtual node per substrate node
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end

    # node capacity : NOT USELESS AHHHHHHHHh
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(v_network[v_node][:dem] * x[v_node, s_node] for v_node in vertices(v_network)) <= sum(s_network[s_node][:cap]))
    end

    ## Edges 
    
    # edge capacity (undirected version !)
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( v_network[src(v_edge), dst(v_edge)][:dem] * (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
                for v_edge in edges(v_network)) 
            <= 
            instance.s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(instance.s_network)
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
    #=
    for s_node in vertices(instance.s_network)
        for v_node in vertices(v_network)
            for v_node_2 in neighbors(v_network, v_node)
                v_edge = get_edge(v_network, v_node, v_node_2)
                if src(v_edge) == v_node # a bit ugly, but this is because the vn is undirected BUT we use it as directed.
                    @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                        >= x[v_node, s_node])
                end
            end
        end
    end
    =#
    # why would this not work ?
    
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end
    
    
    #=
    # Simple path constraints, only useful for porta.
    # Note that non-simple path and subtours are possible with the formulation, 
    # but will never appear in practice due to being expensive for nothing.
    for s_node in vertices(instance.s_network)
        for v_node in vertices(v_network)
            for v_edge in get_out_edges(v_network, v_node)
                @constraint(model, 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node)) 
                    <= 1 - x[v_node, s_node] )
            end
        end
    end
    # to remove loops..
    for v_edge in edges(v_network)
        for s_edge in edges(instance.s_network)
            @constraint(model, y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] 
                + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))] 
                <= 1 )
        end
    end

    =#
    
end



function set_up_problem_placement_restrict(instance, model, placement_restriction)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    ### Variables
    @variable(model, x[vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);

    

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);




    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in placement_restriction[v_node]) == 1)
    end

    # one to one : one virtual node per substrate node
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end

    # node capacity : NOT USELESS AHHHHHHHHh
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(v_network[v_node][:dem] * x[v_node, s_node] for v_node in vertices(v_network)) <= sum(s_network[s_node][:cap]))
    end

    ## Edges 
    
    # edge capacity (undirected version !)
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( v_network[src(v_edge), dst(v_edge)][:dem] * (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
                for v_edge in edges(v_network)) 
            <= 
            instance.s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end

    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end
    
end




function solve_heuristic_relax(instance)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir


    model_prelim = Model(CPLEX.Optimizer)
    set_up_problem(instance, model_prelim)
    relax_integrality(model_prelim)
    print("Starting solving prelims ")
    optimize!(model_prelim)

    x_values = value.(model_prelim[:x])
    placement_restriction = Dict()
    for v_node in vertices(v_network)
        sorted_s_nodes = sort(collect(vertices(s_network_dir)), by=s_node -> -x_values[v_node, s_node])
        placement_restriction[v_node] = sorted_s_nodes[1:30]
    end
    
    # model improved
    model_restrict = Model(CPLEX.Optimizer)
    set_up_problem_placement_restrict(instance, model_restrict, placement_restriction)
    print("Starting solving model_restrict ")
    optimize!(model_restrict)

end





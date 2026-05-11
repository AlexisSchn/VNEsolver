using Revise, JuMP, CPLEX, Gurobi
includet("../utils/import_utils.jl")




# ========== CLASSICAL STUFF

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

    
    # Outgoing edges cap: pretty stupid but useful
    i = 0
    for v_node in vertices(v_network)
        for s_node in vertices(s_network)
            v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
            necessary_bw = 0 + sum(v_network[src(v_edge), dst(v_edge)][:dem] for v_edge in v_edges_incident)

            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = 0 +sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident)
            if necessary_bw > available_bw
                i+=1
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end
    #println("We get this to delete: $i")


    # NON LINEAR STUFF
    # a virtual edge should not be routed both ways bruh
    for v_edge in edges(v_network)
        for s_edge in edges(s_network)
            s_edge_plus = get_edge(s_network_dir, src(s_edge), dst(s_edge))
            s_edge_minus = get_edge(s_network_dir, dst(s_edge), src(s_edge))
            @constraint(model, y[v_edge, s_edge_plus] * y[v_edge, s_edge_minus] == 0)
        end
    end
    
end



function solve_compact(instance, time_solver = 30, stay_silent=false)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir


    model = Model(CPLEX.Optimizer)
    set_up_problem(instance, model)

    #set_optimizer_attribute(model, "CPXPARAM_Emphasis_MIP", 1)

    #relax_integrality(problem.model)
    set_time_limit_sec(model, time_solver)
    if stay_silent
        set_silent(model)
    else
        print("Starting solving... ")
    end
    optimize!(model)

    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("Infeasible or unfinished: $status")
        return
    end

    x_values = value.(model[:x])
    y_values = value.(model[:y])

    if !stay_silent
        println("Node placement:")
        for v_node in vertices(v_network)
            for s_node in vertices(s_network_dir)
                if x_values[v_node, s_node] > 0.5
                    println("$v_node is placed on $s_node")
                end
            end
        end
        println("\nEdge routing:")
        for v_edge in edges(v_network)
            print("Routing of $v_edge : ")
            for s_edge in edges(s_network_dir)
                if y_values[v_edge, s_edge] > 0.5
                    print(" $s_edge")
                end
            end
            print("\n")
        end
    end


end

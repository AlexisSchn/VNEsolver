

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
    
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end
    
    
    
end



function add_star_constraints(instance, model)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    nb_vnodes = length(vertices(v_network))
    # star stuff
    # Get a list of nodes and their degrees
    node_degrees = [(v, degree(v_network, v)) for v in vertices(v_network)]

    # Sort nodes by degree in descending order and take the top five
    dense_v_nodes = sort(node_degrees, by=x -> -x[2])[1:nb_vnodes]




    

    #----------  Local constraints 1
    diff_deg=1
    δ2_edges = Dict()
    for s_node in vertices(s_network_dir)
        edges_in_distance = δ2(s_network_dir, s_node)
        #println("Edge in δ2 of $s_node : $edges_in_distance")
        δ2_edges[s_node] = edges_in_distance
    end

    for (v_node, v_deg) in dense_v_nodes
        v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
        for s_node in vertices(s_network_dir)
            s_deg = degree(instance.s_network, s_node) 
            if s_deg < v_deg
                # Here we need to be extra careful about the pseudo orientation of the virtual network
                trucmuche = @expression(model, 0.)
                for v_edge in v_edges_incident
                    if src(v_edge) == v_node
                        trucmuche += sum(model[:y][v_edge, s_edge]  for s_edge in δ2_edges[s_node])
                    else
                        trucmuche += sum(model[:y][v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  for s_edge in δ2_edges[s_node])
                    end
                end
                if v_deg - s_deg > diff_deg
                    cons = @constraint(model, 
                        trucmuche >= (v_deg - s_deg) * model[:x][v_node, s_node]
                    )
                end
            end
        end
    end



    #----------  Local constraints 2
    diff_deg=1
    for (v_node, v_deg) in dense_v_nodes
        v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
        for s_node in vertices(s_network_dir)
            s_deg = degree(instance.s_network, s_node) 
            if v_deg - s_deg > diff_deg

                # Here we need to be extra careful about the pseudo orientation of the virtual network
                trucmuche = @expression(model, 0.)
                for v_edge in v_edges_incident
                    if src(v_edge) == v_node
                        trucmuche += sum(model[:y][v_edge, s_edge] for neighbor in neighbors(s_network_dir, s_node) for s_edge in get_out_edges(s_network_dir, neighbor) )
                    else
                        trucmuche += sum(model[:y][v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  for neighbor in neighbors(s_network_dir, s_node) for s_edge in get_out_edges(s_network_dir, neighbor) )
                    end
                end
                
                @constraint(model, 
                    trucmuche >= (v_deg - s_deg) * model[:x][v_node, s_node] + sum(v_deg * model[:x][v_node, s_neighbor] for s_neighbor in neighbors(s_network_dir, s_node))
                )
            end
        end
    end




    #---------------- global constraints + precise
    for (v_node, v_deg) in dense_v_nodes
        trucmuche = @expression(model, v_deg)
        for s_node in vertices(s_network_dir)
            trucmuche += sum(model[:x][v_node, s_node] * (v_deg - degree(instance.s_network, s_node)))
        end

        v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
        
        @constraint(model, 
            trucmuche <= sum(model[:y][v_edge, s_edge] for v_edge in v_edges_incident for s_edge in edges(s_network_dir))
        )
    end


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
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1) ### A RAJOUTER ABSOLUMNET ! (ou alors tu ne créés pas les autres variables).
        # sinon il est pas au courant qu'il peut "supprimer" les x pas valables...
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



function solve_compact_pl_1vn(instance, time_solver = 30)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir


    model_prelim = Model(CPLEX.Optimizer)
    set_up_problem(instance, model_prelim)
    add_star_constraints(instance, model_prelim)


    relax_integrality(model_prelim)
    print("Starting solving prelims ")
    optimize!(model_prelim)


    x_values = value.(model_prelim[:x])
    y_values = value.(model_prelim[:y])

    println("Node placement:")
    for v_node in vertices(v_network)
        println("   Node $v_node:")
        for s_node in vertices(s_network_dir)
            if x_values[v_node, s_node] > 0.001
                println("       $s_node : $(x_values[v_node, s_node])")
            end
        end
    end

    println("\nEdge routing:")
    for v_edge in edges(v_network)
        println("   $v_edge : ")
        for s_edge in edges(s_network_dir)
            if y_values[ v_edge, s_edge] > 0.001
                println("         $s_edge : $(y_values[v_edge, s_edge])")
            end
        end
        print("\n")
    end


end




function solve_heuristic_relax(instance)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir


    model_prelim = Model(CPLEX.Optimizer)
    set_up_problem(instance, model_prelim)
    add_star_constraints(instance, model_prelim)


    relax_integrality(model_prelim)
    print("Starting solving prelims ")
    optimize!(model_prelim)

    x_values = value.(model_prelim[:x])
    placement_restriction = Dict()
    for v_node in vertices(v_network)
        sorted_s_nodes = sort(collect(vertices(s_network_dir)), by=s_node -> -x_values[v_node, s_node])
        placement_restriction[v_node] = sorted_s_nodes[1:12]
    end
    
    # model improved
    model_restrict = Model(CPLEX.Optimizer)
    set_up_problem_placement_restrict(instance, model_restrict, placement_restriction)
    print("Starting solving model_restrict ")
    optimize!(model_restrict)

end






# edge in delta^2(u)
function δ2(g, node)
    δ2_edges = []
    neighbors_node = neighbors(g, node)
    for neighbor in neighbors_node
        #push!(edges_in_distance, get_edge(g, origin_node, neighbor))
        #push!(edges_in_distance, get_edge(g, neighbor, origin_node))
        neighbors_of_neighbor = neighbors(g, neighbor)
        for neighbor_of_neighbor in neighbors_of_neighbor
            if (neighbor_of_neighbor != node) && !(neighbor_of_neighbor in neighbors_node)
                push!(δ2_edges, get_edge(g, neighbor, neighbor_of_neighbor))
                #push!(edges_in_distance, get_edge(g, neighbor_of_neighbor, neighbor))
            end
        end
    end
    return δ2_edges
end
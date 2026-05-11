includet("compact_undir.jl")


function solve_compact_pretreatment_star(instance, is_integer = true, time_solver = 30)

        
    model = Model(CPLEX.Optimizer)
    set_up_problem_undir_1vn_1t1(instance, model)
    set_time_limit_sec(model, time_solver)
    #optimize!(problem.model)



    # star stuff
    # Get a list of nodes and their degrees
    v_network = instance.v_network
    node_degrees = [(v, degree(v_network, v)) for v in vertices(v_network)]

    # Sort nodes by degree in descending order and take the top five
    nb_v_node = 2
    top_five_nodes = sort(node_degrees, by=x -> -x[2])[1:nb_v_node]


    close_edges = []
    for s_node in vertices(s_network_dir)
        edges_in_distance = edges_within_hop(s_network_dir, s_node)
        #println("Edge close to $s_node : $edges_in_distance")
        push!(close_edges, edges_in_distance)
    end

    nb_cons = 0
    for (v_node, v_deg) in top_five_nodes
        v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
        for s_node in vertices(s_network_dir)
            s_deg = degree(instance.s_network, s_node) 
            if s_deg < v_deg

                # Here we need to be extra careful about the pseudo orientation of the virtual network
                trucmuche = @expression(problem.model, 0.)
                for v_edge in v_edges_incident
                    if src(v_edge) == v_node
                        trucmuche += sum(problem.model[:y][v_edge, s_edge]  for s_edge in close_edges[s_node])
                    else
                        trucmuche += sum(problem.model[:y][v_network, v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  for s_edge in close_edges[s_node])
                    end
                end
                cons = @constraint(problem.model, 
                    trucmuche >= (2 * v_deg - s_deg) * problem.model[:x][v_network, v_node, s_node]
                )

                nb_cons += 1
                #println("Wow look at this constraint: $cons")

            end
        end
    end


    println("\nWe have added $nb_cons star cuts")
    is_integer = false  
    if !is_integer
        relax_integrality(model)
    end
    optimize!(model)

    x_values = value.(model[:x])
    y_values = value.(model[:y])


    if is_integer
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
    else
        v_network = instance.v_networks[1]
        println("Node placement:")
        for v_node in vertices(v_network)
            println("   Node $v_node:")
            for s_node in vertices(s_network_dir)
                if x_values[v_network, v_node, s_node] > 0.001
                    println("       $s_node : $(x_values[v_network, v_node, s_node])")
                end
            end
        end
        println("\nEdge routing:")
        for v_edge in edges(v_network)
            println("   $v_edge : ")
            for s_edge in edges(s_network_dir)
                if y_values[v_network, v_edge, s_edge] > 0.001
                    println("         $s_edge : $(y_values[v_network, v_edge, s_edge])")
                end
            end
            print("\n")
        end
    end


end


function edges_within_hop(g, origin_node)
    edges_in_distance = []
    neighbors_node = neighbors(g, origin_node)
    for neighbor in neighbors_node
        push!(edges_in_distance, get_edge(g, origin_node, neighbor))
        #push!(edges_in_distance, get_edge(g, neighbor, origin_node))

        neighbors_of_neighbor = neighbors(g, neighbor)
        for neighbor_of_neighbor in neighbors_of_neighbor
            if neighbor_of_neighbor != origin_node
                push!(edges_in_distance, get_edge(g, neighbor, neighbor_of_neighbor))
                #push!(edges_in_distance, get_edge(g, neighbor_of_neighbor, neighbor))
            end
        end
    end
    return edges_in_distance
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


function edges_within_two_hops(g::Graph, node)
    # Collect all edges within distance 2 from the given node
    edges_in_distance = Set{Edge}()

    # First-level neighbors (distance 1)
    neighbors1 = neighbors(g, node)
    for neighbor1 in neighbors1
        edges_in_distance ∪ Set(edges_incident(g, neighbor1))  # Add edges to/from first neighbors
        
        # Second-level neighbors (distance 2)
        neighbors2 = neighbors(g, neighbor1)
        for neighbor2 in neighbors2
            edges_in_distance ∪ Set(edges_incident(g, neighbor2))  # Add edges to/from second neighbors
        end
    end

    return edges_in_distance
end



function solve_compact_pretreatment_star_v2(instance, is_integer = true, nb_nodes = 5, time_solver = 30)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
        
    model = Model(CPLEX.Optimizer)
    set_up_problem_undir_1vn_1t1(instance, model)
    set_time_limit_sec(model, time_solver)
    #optimize!(problem.model)



    # star stuff
    # Get a list of nodes and their degrees
    node_degrees = [(v, degree(v_network, v)) for v in vertices(v_network)]

    # Sort nodes by degree in descending order and take the top five
    dense_v_nodes = sort(node_degrees, by=x -> -x[2])[1:nb_nodes]



    #---------- 1: Local constraints
    δ2_edges = Dict()
    for s_node in vertices(s_network_dir)
        edges_in_distance = δ2(s_network_dir, s_node)
        #println("Edge in δ2 of $s_node : $edges_in_distance")
        δ2_edges[s_node] = edges_in_distance
    end

    nb_cons = 0
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
                add_locals = true      
                if add_locals
                    if v_deg - s_deg > 1
                        cons = @constraint(model, 
                            trucmuche >= (v_deg - s_deg) * model[:x][v_node, s_node]
                        )
                    end
                end
                
                
                nb_cons += 1
            end
        end
    end
    println("\nWe have added $nb_cons local star cuts\n\n\n")



    #---------------- 2 : global constraints
    nb_cons_2 = 0
    for (v_node, v_deg) in dense_v_nodes

        trucmuche = @expression(model, v_deg)
        for s_node in vertices(s_network_dir)
            trucmuche += sum(model[:x][v_node, s_node] * (v_deg - degree(instance.s_network, s_node)))
        end


        v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
        
        add_globals = true 
        if add_globals
            cons = @constraint(model, 
                trucmuche <= sum(model[:y][v_edge, s_edge] for v_edge in v_edges_incident for s_edge in edges(s_network_dir))
            )
        end
        
        
        #println("New cons: $cons")
        nb_cons_2 += 1
    end
    println("\nWe have added $(nb_cons_2) global star cuts\n\n\n")

    #---------- solving
    if !is_integer
        relax_integrality(model)
    end
    optimize!(model)

    println("\n\nWell the objective value is : $(objective_value(model))\n\n")

    x_values = value.(model[:x])
    y_values = value.(model[:y])

    if is_integer
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
    else
        println("Node placement:")
        for v_node in vertices(v_network)
            println("   Node $v_node:")
            for s_node in vertices(s_network_dir)
                if x_values[v_node, s_node] > 0.001
                    #println("       $s_node : $(x_values[v_node, s_node])")
                end
            end
        end

        println("\nEdge routing: too long")
        #=
        for v_edge in edges(v_network)
            println("   $v_edge : ")
            for s_edge in edges(s_network_dir)
                if y_values[v_edge, s_edge] > 0.001
                    println("         $s_edge : $(y_values[v_edge, s_edge])")
                end
            end
            print("\n")
        end
        =#
    end


end




function solve_compact_basic(instance, is_integer = true, time_solver = 30)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
        
    model = Model(CPLEX.Optimizer)
    set_up_problem_undir_1vn_1t1(instance, model)
    set_time_limit_sec(model, time_solver)
    #optimize!(problem.model)





    #---------- solving
    if !is_integer
        relax_integrality(model)
    end
    optimize!(model)

    println("\n\nWell the objective value is : $(objective_value(model))\n\n")

    x_values = value.(model[:x])
    y_values = value.(model[:y])

    if is_integer
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
    else
        println("Node placement:")
        for v_node in vertices(v_network)
            println("   Node $v_node:")
            for s_node in vertices(s_network_dir)
                if x_values[v_node, s_node] > 0.001
                    println("       $s_node : $(x_values[v_node, s_node])")
                end
            end
        end

        println("\nEdge routing: too long")
        #=
        for v_edge in edges(v_network)
            println("   $v_edge : ")
            for s_edge in edges(s_network_dir)
                if y_values[v_edge, s_edge] > 0.001
                    println("         $s_edge : $(y_values[v_edge, s_edge])")
                end
            end
            print("\n")
        end
        =#
    end


end


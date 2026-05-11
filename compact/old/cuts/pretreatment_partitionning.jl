includet("compact_undir.jl")




# Ne marche qu'en one to one pour l'instant...
function solve_compact_pretreatment_partitionning(instance, node_partitionning, time_solver = 30)
    

    # We have to create a digraph with symetric edges yaknow
    s_network_dir = instance.s_network_dir
    v_network = instance.v_network
    s_network = instance.s_network
        
    model = Model(CPLEX.Optimizer)

    #set_optimizer_attribute(model, "CPXPARAM_MIP_Strategy_VariableSelect", 4)
    set_up_problem_undir_1vn_1t1(instance, model)

    # Subgraph stuff....
    for (i_subgraph, v_nodes) in enumerate(node_partitionning)
        print("For subgraph $i_subgraph... ")
        subgraph = my_induced_subgraph(v_network, v_nodes, "subgraph_$i_subgraph")
        print_graph(subgraph)
        # Résolution
        subinstance = Instance_Undir_VNE_1s(subgraph, s_network)
        submodel = Model(CPLEX.Optimizer)
        set_silent(submodel)
        set_up_problem_undir_1vn_1t1(subinstance, submodel)
        set_time_limit_sec(submodel, time_solver)
        optimize!(submodel)
        min_routing = objective_value(submodel)

        print("Min edges: $min_routing instead of $(length(edges(subgraph)))\n")

        #adding the constraints
        @constraint(model,
            sum( model[:y][get_edge(v_network, v_nodes[src(v_edge)], v_nodes[dst(v_edge)]), s_edge] for v_edge in edges(subgraph) for s_edge in edges(s_network_dir))
            >=
            min_routing)
            
    end 

    

    # Solving
    print("Starting solving... ")
    set_time_limit_sec(model, time_solver)
    
    is_integer = false
    if !is_integer
        relax_integrality(model)
    end

    optimize!(model)
    println("Finale value: $(objective_value(model))")

    return

    x_values = value.(problem.model[:x])
    y_values = value.(problem.model[:y])

    v_network = instance.v_networks[1]

    if is_integer
        println("Node placement:")
        for v_node in vertices(v_network)
            for s_node in vertices(s_network_dir)
                if x_values[v_network, v_node, s_node] > 0.5
                    println("$v_node is placed on $s_node")
                end
            end
        end
        println("\nEdge routing:")
        for v_edge in edges(v_network)
            print("Routing of $v_edge : ")
            for s_edge in edges(s_network_dir)
                if y_values[v_network, v_edge, s_edge] > 0.5
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



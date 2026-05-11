includet("../../../utils/import_utils.jl")
includet("solution_path_fractional.jl")

using JuMP, CPLEX, Random


struct Master_Problem
    model
    x_variables
    lambda_variables
end


function number_of_column(master_problem)

    num_cols = 0
    
    for v_n in keys(master_problem.lambda_variables)
        for v_edge in keys(master_problem.lambda_variables[v_n])
            num_cols += length(master_problem.lambda_variables[v_n][v_edge])
        end
    end
    
    return num_cols 
end


function set_up_master_problem(instance)

    print("Constructing empty model...")
    #### Model
    model = Model(CPLEX.Optimizer)
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)
    
    # Variables
    x_variables = @variable(model, x[v_network in instance.v_networks, vertices(v_network), vertices(instance.s_network)], binary=true)

    lambda_variables = Dict()
    for v_network in instance.v_networks
        lambda_variables[v_network] = Dict()
        for v_edge in edges(v_network)
            lambda_variables[v_network][v_edge] = Dict()
        end
    end


    # Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_network, v_node, s_node] 
        for v_network in instance.v_networks for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    @objective(model, Min, placement_cost );



    # Constraints

    #### NODES

    # one substrate node per virtual node
    for v_network in instance.v_networks
        for v_node in vertices(v_network)
            @constraint(model, sum(x[v_network, v_node, :]) == 1)
        end
    end

    
    # if one to one : one virtual node per substrate node
    one_to_one = true
    if one_to_one 
        for s_node in vertices(instance.s_network)
            for v_network in instance.v_networks
                @constraint(model, sum(x[v_network, v_node, s_node] for v_node in vertices(v_network)) <= 1)
            end
        end
    end

    #capacity
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(  v_network[v_node][:dem] * x[v_network, v_node, s_node] 
                                for v_network in instance.v_networks for v_node in vertices(v_network) ) 
                    <= instance.s_network[s_node][:cap] )
    end
    
    
    #### EDGES

    # one path per v_edge
    @constraint(
        model, 
        path_selec[v_network in instance.v_networks, v_edge in edges(v_network)],
        0 == 1
    );

    # capacity
    @constraint(
        model,
        capacity_s_edge[s_edge in edges(instance.s_network)],
        0 <= instance.s_network[src(s_edge), dst(s_edge)][:cap]  
    );

    # start
    @constraint(
        model, 
        start[v_network in instance.v_networks, v_edge in edges(v_network), s_node in vertices(instance.s_network)],
        0 == x[v_network, src(v_edge), s_node]
    );
        
    # terminus
    @constraint(
        model, 
        destination[v_network in instance.v_networks, v_edge in edges(v_network), s_node in vertices(instance.s_network)],
        0 == x[v_network, dst(v_edge), s_node]
    );



    println(" done.")

    # relax_integrality(model);
    return(Master_Problem(model, x_variables, lambda_variables))
end



struct PricerProblem
    model
    s_network
    v_network
    v_edge
end

function set_up_pricer(s_network, v_network, v_edge)
    #### Model
    model = Model(CPLEX.Optimizer)
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)

    ### Variables
    @variable(model, start[vertices(s_network)], binary=true);
    @variable(model, terminus[vertices(s_network)], binary=true);
    @variable(model, y[edges(s_network)], binary=true);

    @objective(model, Min, 0); # willl be updated throughout resolution

    ### Constraints

    ## Nodes

    # one start and one terminus
    @constraint(model, sum(start[:]) == 1)
    @constraint(model, sum(terminus[:]) == 1)

    # capacity
    for s_node in vertices(s_network)
        @constraint(model, v_network[src(v_edge)][:dem] * start[s_node] + v_network[dst(v_edge)][:dem] * terminus[s_node] <= s_network[s_node][:cap])
    end

    one_to_one = true
    if one_to_one
        @constraint(model, start .+ terminus .<= 1)
    end


    ## Edges

    # flow constraints
    for s_node in vertices(s_network)
        @constraint(model, 
            start[s_node] + sum(y[s_edge] for s_edge in get_in_edges(s_network, s_node)) 
            == 
            sum(y[s_edge] for s_edge in get_out_edges(s_network, s_node)) + terminus[s_node] )
    end

    # capacity
    for s_edge in edges(s_network)
        if s_network[src(s_edge), dst(s_edge)][:cap] < v_network[src(v_edge), dst(v_edge)][:dem]
            @constraint(model, y[s_edge] == 0)
        end
    end

    return PricerProblem(model, s_network, v_network, v_edge)
end


function update_and_solve_pricer(pricer, pi, beta, alphas, alphat)
    
    
    #updating ojective values with dual variables
    @objective(pricer.model, Min, -pi[pricer.v_network][pricer.v_edge])

    for s_node in vertices(pricer.s_network)
        set_objective_coefficient(pricer.model, pricer.model[:start][s_node], - alphas[pricer.v_network][pricer.v_edge][s_node])
        set_objective_coefficient(pricer.model, pricer.model[:terminus][s_node], - alphat[pricer.v_network][pricer.v_edge][s_node])
    end

    for s_edge in edges(pricer.s_network)
        set_objective_coefficient(pricer.model, pricer.model[:y][s_edge], (pricer.s_network[src(s_edge), dst(s_edge)][:cost] - beta[s_edge]) * 
                            pricer.v_network[src(pricer.v_edge), dst(pricer.v_edge)][:dem])
    end
    set_silent(pricer.model)
    optimize!(pricer.model)

    # Get the solution
    start_values = value.(pricer.model[:start])
    terminus_values = value.(pricer.model[:terminus])
    y_values = value.(pricer.model[:y])

    start_node = 0
    terminus_node = 0
    for s_node in vertices(pricer.s_network)
        if start_values[s_node] > 0.99
            start_node = s_node
        end
        if terminus_values[s_node] > 0.99
            terminus_node = s_node
        end
    end

    used_edges = []
    for s_edge in edges(pricer.s_network)
        if y_values[s_edge] > 0.99
            push!(used_edges, s_edge)
        end
    end
    path = order_path(pricer.s_network, used_edges, start_node, terminus_node) 

    return path, objective_value(pricer.model)
end

function add_path(master_problem, v_network, v_edge, s_path)
    variable_path = @variable(master_problem.model, binary = true)
    master_problem.lambda_variables[v_network][v_edge][s_path] = variable_path
    set_objective_coefficient(master_problem.model, variable_path, v_network[src(v_edge), dst(v_edge)][:dem] * s_path.cost)
    set_normalized_coefficient(master_problem.model[:path_selec][v_network, v_edge], variable_path, 1)
    for s_edge in s_path.edges
        set_normalized_coefficient(master_problem.model[:capacity_s_edge][s_edge], variable_path, v_network[src(v_edge), dst(v_edge)][:dem])
    end
    set_normalized_coefficient(master_problem.model[:start][v_network, v_edge, s_path.src], variable_path, 1)
    set_normalized_coefficient(master_problem.model[:destination][v_network, v_edge, s_path.dst], variable_path, 1)  
end



function get_duals(master_problem)
    beta = Dict()
    for s_edge in edges(instance.s_network)
        beta[s_edge] =  dual(master_problem.model[:capacity_s_edge][s_edge])
    end
    
    pi = Dict()
    for v_network in instance.v_networks
        pi[v_network] = Dict()
        for v_edge in edges(v_network)
            pi[v_network][v_edge] = dual(master_problem.model[:path_selec][v_network, v_edge])
        end
    end
    
    alphas = Dict()
    for v_network in instance.v_networks
        alphas[v_network] = Dict()
        for v_edge in edges(v_network)
            alphas[v_network][v_edge] = []
            for s_node in vertices(instance.s_network)
                push!(alphas[v_network][v_edge], dual(master_problem.model[:start][v_network, v_edge, s_node]))
            end
        end
    end
    
    alphat = Dict()
    for v_network in instance.v_networks
        alphat[v_network] = Dict()
        for v_edge in edges(v_network)
            alphat[v_network][v_edge] = []
            for s_node in vertices(instance.s_network)
                push!(alphat[v_network][v_edge], dual(master_problem.model[:destination][v_network, v_edge, s_node]))
            end
        end
    end    

    return pi, beta, alphas, alphat
end




function solve_master_problem(master_problem)
    set_silent(master_problem.model)
    relax_integrality(master_problem.model)
    optimize!(master_problem.model)
    println("Valeur objective : " * string(objective_value(master_problem.model)))
end

function solve_instance_path_int(instance)
    # set up master problem and satelites problem
    master_problem = set_up_master_problem(instance)


    # generate paths
    print("Generating all paths...")
    paths_substrate = get_shortest_paths(instance.s_network, 1)
    all_paths = []
    for pair_nodes in keys(paths_substrate)
        for path in paths_substrate[pair_nodes]
            push!(all_paths, path)
        end
    end
    println(" done")
    
    for v_network in instance.v_networks
        for v_edge in edges(v_network)
            for s_path in all_paths
                add_path(master_problem, v_network, v_edge, s_path)
            end
        end
    end

    optimize!(master_problem.model)

    println("Solving done ! Number of columns : " * string(number_of_column(master_problem)))

    x_values = value.(master_problem.x_variables)
    mappings = MappingClassic[]
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            for s_node in vertices(instance.s_network)
                if x_values[v_network, v_node, s_node] >= 0.9
                    append!(node_placement, s_node)
                end
            end
        end
        edge_routing = Dict()
        for v_edge in edges(v_network)
            for s_path in keys(master_problem.lambda_variables[v_network][v_edge])
                if value.(master_problem.lambda_variables[v_network][v_edge][s_path]) >= 0.9
                    edge_routing[v_edge] = s_path
                end
            end
        end

        m = MappingClassic(v_network, instance.s_network, node_placement, edge_routing)
        push!(mappings, m)
    end

    return mappings
end    

function solve_relax_with_all_columns(instance)

    # set up master problem and satelites problem
    master_problem = set_up_master_problem(instance)


    # generate initial set of columns
    choice_gene_path = 2
    if choice_gene_path == 1
        for s_node in vertices(instance.s_network)
            s_path = Path(s_node, s_node, [], 0)
            for v_network in instance.v_networks
                for v_edge in edges(v_network)
                    add_path(master_problem, v_network, v_edge, s_path)
                end
            end
        end
    else
        print("Generating all paths...")
        paths_substrate = get_shortest_paths(instance.s_network, 1)
        all_paths = []
        for pair_nodes in keys(paths_substrate)
            for path in paths_substrate[pair_nodes]
                
                one_to_one = false
                if path.src != path.dst && one_to_one
                    push!(all_paths, path)
                elseif !one_to_one
                    push!(all_paths, path)
                end

            end
        end
        println(" done")
        
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                for s_path in all_paths
                    add_path(master_problem, v_network, v_edge, s_path)
                end
            end
        end
    
    end

    solve_master_problem(master_problem)

    println("Solving done ! Number of columns : " * string(number_of_column(master_problem)))


    # getting the fractional solution
    x_values = value.(master_problem.x_variables)
    mappings = Mapping_Path_Fractional[]
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            push!(node_placement, [])
            for s_node in vertices(instance.s_network)
                push!(node_placement[v_node], x_values[v_network, v_node, s_node])
            end
        end
        edge_routing = Dict()
        for v_edge in edges(v_network)
            edge_routing[v_edge] = Dict()
            for s_path in keys(master_problem.lambda_variables[v_network][v_edge])
                edge_routing[v_edge][s_path] = value.(master_problem.lambda_variables[v_network][v_edge][s_path])
            end
        end

        m = Mapping_Path_Fractional(v_network, instance.s_network, node_placement, edge_routing)
        push!(mappings, m)
    end

    return mappings

end


function solve_column_generation(instance)

    # set up master problem and pricer problem
    master_problem = set_up_master_problem(instance)
    pricers_problem = Dict()
    for v_network in instance.v_networks
        pricers_problem[v_network] = Dict()
        for v_edge in edges(v_network)
            pricers_problem[v_network][v_edge] = set_up_pricer(instance.s_network, v_network, v_edge)
        end
    end

    # generate initial set of columns
    print("Generating all paths...")
    paths_substrate = get_shortest_paths(instance.s_network, 1)
    all_paths = []
    paths_of_v_edge = Dict()
    for v_network in instance.v_networks
        paths_of_v_edge[v_network] = Dict()
        for v_edge in edges(v_network)
            paths_of_v_edge[v_network][v_edge] = []
        end
    end
    #=
    Random.seed!(1)
    for pair_nodes in keys(paths_substrate)
        for path in paths_substrate[pair_nodes]
            
            one_to_one = true
            if path.src != path.dst && one_to_one
                
                i = rand((1, 10))
                imlucky = false
                if i > 9
                    imlucky = true
                end
                
                if imlucky
                    push!(all_paths, path)
                end
            elseif !one_to_one
                push!(all_paths, path)
            end

        end
    end
    =#
    Random.seed!(1)
    for pair_nodes in keys(paths_substrate)
        for path in paths_substrate[pair_nodes]
            
            if path.src != path.dst 
                push!(all_paths, path)
            end

        end
    end
    println(" done")
    
    for v_network in instance.v_networks
        for v_edge in edges(v_network)
            for s_path in all_paths
                add_path(master_problem, v_network, v_edge, s_path)
            end
        end
    end
    
    keep_on_trying = true
    while keep_on_trying
        keep_on_trying = false
        solve_master_problem(master_problem)

        pi, beta, alphas, alphat = get_duals(master_problem)

        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                path, value = update_and_solve_pricer(pricers_problem[v_network][v_edge], pi, beta, alphas, alphat)
                if value < 0
                    if !path_in_paths(path, paths_of_v_edge[v_network][v_edge])
                        push!(paths_of_v_edge[v_network][v_edge], path)
                        add_path(master_problem, v_network, v_edge, path)
                        #println("new path added for edge "* string(v_edge) * " : " * string(path))
                        keep_on_trying = true
                    else 
                        #println("Path already existed... value: " * string(value))
                    end
                end
            end
        end
    end

    println("Column generation over. Number of columns : " * string(number_of_column(master_problem)))
    

    # getting the fractional solution
    x_values = value.(master_problem.x_variables)
    mappings = Mapping_Path_Fractional[]
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            push!(node_placement, [])
            for s_node in vertices(instance.s_network)
                push!(node_placement[v_node], x_values[v_network, v_node, s_node])
            end
        end
        edge_routing = Dict()
        for v_edge in edges(v_network)
            edge_routing[v_edge] = Dict()
            for s_path in keys(master_problem.lambda_variables[v_network][v_edge])
                edge_routing[v_edge][s_path] = value.(master_problem.lambda_variables[v_network][v_edge][s_path])
            end
        end

        m = Mapping_Path_Fractional(v_network, instance.s_network, node_placement, edge_routing)
        push!(mappings, m)
    end

    return mappings

end




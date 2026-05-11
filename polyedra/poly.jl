
using JuMP, CPLEX, Gurobi
using Polyhedra, CDDLib, XPORTA


includet("../utils/import_utils.jl")



function set_up_problem(instance, model)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    ### Variables
    @variable(model, x[vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);

    

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * x[v_node, s_node] 
        for v_node in vertices(v_network), s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * y[v_edge, s_edge]
        for v_edge in edges(v_network), s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);




    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    @constraint(model, [v_node in vertices(v_network)], 
        sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1
    )

    # one to one : one virtual node per substrate node
    @constraint(model, [s_node in vertices(s_network)], 
        sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1
    )
    
    # node capacity : NOT USELESS AHHHHHHHHh
    @constraint(model, [s_node in vertices(s_network)], 
        sum( x[v_node, s_node] for v_node in vertices(v_network)) <= sum(s_network[s_node][:cap])
    )


    ## Edges 
    
    # edge capacity (undirected version !)
    @constraint(model, [s_edge in edges(instance.s_network)],
        sum( y[v_edge, s_edge] + y[v_edge, get_reverse_edge(s_network_dir, s_edge)]  for v_edge in edges(v_network)) 
            <= instance.s_network[src(s_edge), dst(s_edge)][:cap] 
    )
    
    # Flow conservation
    @constraint(model, [ s_node in vertices(instance.s_network), v_edge in edges(v_network)],
        x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
        ==
        sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
            sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
    )

    
    #= Departure constraints - useless for porta
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end
    =#
    
    
    # TODO: just do not take in account the solution like this? Cuze you need subtour constraint, which are annoying af.
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
    # to remove loops.. but here, your not removing all of them.
    for v_edge in edges(v_network)
        for s_edge in edges(instance.s_network)
            @constraint(model, y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] 
                + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))] 
                <= 1 )
        end
    end
    
    
    
end



# Get all solutions of the instance
function get_sols(instance)

    model = Model(Gurobi.Optimizer)
    set_silent(model)
    s_network_dir = instance.s_network_dir

    set_up_problem(instance, model)
    println("Presolve at 0 - potential issue?")
    set_optimizer_attribute(model, "PoolSearchMode", 2)
    set_optimizer_attribute(model, "PoolSolutions", 10000)
    
    print("Model set, starting to solve... ")

    optimize!(model)
    #solution_summary(model)
    nb_solutions = result_count(model)
    println("Done, there are $nb_solutions solutions")
    sols = Vector{Vector{Int64}}()
    for i in 1:nb_solutions
        cur_sol = Vector{Int64}()
        for v_node in vertices(instance.v_network), s_node in vertices(s_network_dir)
            push!(cur_sol, round.(value.(model[:x][v_node, s_node]; result = i)))     
        end

        for v_edge in edges(instance.v_network), s_edge in edges(s_network_dir)
            push!(cur_sol, round.(value.(model[:y][v_edge, s_edge]; result = i)))     
        end

        push!(sols, cur_sol)
    end

    name_variables = []
    for v_node in vertices(instance.v_network), s_node in vertices(instance.s_network)
        push!(name_variables, "x_" * string(v_node) * "_" * string(s_node))
    end
    
    for v_edge in edges(instance.v_network), s_edge in edges(s_network_dir)
        push!(name_variables, "y_" * string(src(v_edge)) * string(dst(v_edge)) * "_" * string(src(s_edge)) * string(dst(s_edge)))
    end

    return sols, name_variables
end


# Get h representation from v representation
function get_hrep(sols, names)
    v_rep = vrep(sols);
    poly = polyhedron(v_rep, CDDLib.Library());
    removehredundancy!(poly)
    println("There is $(length(sols)) solutions and $(length(names)) variables in your instance.")
    println("Computing the hrep... This might take some time... ")
    h_rep = hrep(poly);
    print("finished.")
    return h_rep
end



# Get the dominant - didnt use it really, tried it and it gave weird results...?
function get_dominant_hrep(sols, names)
    nb_variables = length(names)
    rays = Vector{Vector{Int}}()
    for i_var in 1:nb_variables
        ray = zeros(Int64, nb_variables)
        ray[i_var] = 1
        push!(rays, ray)
    end
    v_rep = convexhull(sols...) + conichull(rays...) # This is the "minkovsky sum" or something like that
    poly = polyhedron(v_rep, CDDLib.Library());
    #print(poly)
    removehredundancy!(poly)
    println("There is $(length(sols)) solutions and $(length(names)) variables in your instance.")
    println("Computing the hrep... This might take some time... ")
    h_rep = hrep(poly);

    return(h_rep)

end




# ------ Priting functions


function print_solutions(sols, names)

    for (i_sol, sol) in enumerate(sols)
        println("Solution num $(i_sol)")
        for (i_var, value) in enumerate(sol) 
            if value > 0.5
                println("$(names[i_var])")
            end
        end
    end
end


function print_polytope(hr, names_variables)
    println("There are " * string(length(names_variables)) * " variables.\n")
    println("There are " * string(length(hyperplanes(hr))) * " hyperplanes")

    for h in hyperplanes(hr)
        str_right = ""
        str_left = ""
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                str_left *= "+ $(floor(Int, h.a[i_var]))$(names_variables[i_var])  "
            elseif (h.a[i_var] < -0.0001)
                str_right *= "+ $(floor(Int, -h.a[i_var]))$(names_variables[i_var])  "
            end
        end

        if  h.β > 0.00001
            str_right = " $(floor(Int,h.β)) " * str_right
        elseif h.β < - 0.00001
            str_left = " $(floor(Int,-h.β)) " * str_left
        end
        
        if length(str_left) == 0
            str_left = " 0 "
        end
        if length(str_right) == 0
            str_right = " 0 "
        end

        print(str_left)
        print(" = ")
        println(str_right)
    end

    println("\n\n There are " * string(length(halfspaces(hr))) * " halfspaces")
    i_hs = 1
    for h in halfspaces(hr)
        str_right = ""
        str_left = ""
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                str_left *= "+ $(floor(Int, h.a[i_var]))$(names_variables[i_var])  "
            elseif (h.a[i_var] < -0.0001)
                str_right *= "+ $(floor(Int, -h.a[i_var]))$(names_variables[i_var])  "
            end
        end

        if  h.β > 0.00001
            str_right = " $(floor(Int,h.β)) " * str_right
        elseif h.β < - 0.00001
            str_left = " $(floor(Int,-h.β)) " * str_left
        end
        
        if length(str_left) == 0
            str_left = " 0 "
        end
        if length(str_right) == 0
            str_right = " 0 "
        end
        print(str_right)
        print(" ≥ ")
        print(str_left)
        println("\t($(i_hs))")
        i_hs += 1
    end

end




# This was to write all the constraint to add them to the model directly, by hand, to try and check which one would improve the bound.
function print_constraints_jump(hr, instance)

    for h in halfspaces(hr)
        i_var = 1

        print("@constraint(model, ")
        
        for v_node in vertices(instance.v_network)
            for s_node in vertices(instance.s_network)
                if abs((h.a[i_var])) > 0.01
                    print("+ $(floor(Int, h.a[i_var])) * model[:x][$v_node, $s_node] ")
                end
                i_var += 1
            end
        end
        for v_edge in edges(instance.v_network)
            for s_edge in edges(instance.s_network_dir)
                if (abs(h.a[i_var]) > 0.01)
                    print("+ $(floor(Int, h.a[i_var])) * model[:y][get_edge(instance.v_network, $(src(v_edge)), $(dst(v_edge))), get_edge(instance.s_network_dir, $(src(s_edge)), $(dst(s_edge)))] ")
                end
                i_var += 1
            end
        end
        print(" <= +")
        print(floor(Int,h.β))
        println(")")


    end

end


# Given a solution, would check which facets are violated. But, in my case, since I have hyperplanes, it didnt worked so I abandonned it.
function print_violated_constraints(hr, names, solution) 

    println("\n\n There are " * string(length(halfspaces(hr))) * " halfspaces")

    for h in halfspaces(hr)
        sum_vals = 0
        for i_var in 1:length(names)
            sum_vals += h.a[i_var] * solution[i_var]
        end

        if sum_vals > h.β
            println("New constraint violated ! $sum_vals > $(h.β)")
            print_constraint(h, names)
        end
        
        #println("Ok lets check...")
        #println("Here is my solution: $solution")
        #print("Here is my constraint:  ")
        #print_constraint(h, names)
        #println("So lhs: $(sum_vals), and rhs: $(h.β)")
        
    end

    println("Finished...")

end


function print_constraint_better(h, names_variables) 
    
    str_right = ""
    str_left = ""
    for i_var in 1:length(names_variables)
        if (h.a[i_var] > 0.0001)
            str_left *= "+ $(floor(Int, h.a[i_var]))$(names_variables[i_var])  "
        elseif (h.a[i_var] < -0.0001)
            str_right *= "+ $(floor(Int, -h.a[i_var]))$(names_variables[i_var])  "
        end
    end

    if  h.β > 0.00001
        str_right = " $(floor(Int,h.β)) " * str_right
    elseif h.β < - 0.00001
        str_left = " $(floor(Int,-h.β)) " * str_left
    end
    
    if length(str_left) == 0
        str_left = " 0 "
    end
    if length(str_right) == 0
        str_right = " 0 "
    end
    print(str_right)
    print(" ≥ ")
    println(str_left)


end


function print_constraint(h, names_variables) 
    
    for i_var in 1:length(names_variables)
        if (h.a[i_var] > 0.0001)
            print("+ ")
            print(h.a[i_var])
            print(" ")
            print(names_variables[i_var])
            print("\t")
        elseif (h.a[i_var] < -0.0001)
            print("- ")
            print(-h.a[i_var])
            print(" ")
            print(names_variables[i_var])
            print("\t")
        end
    end
    print(" ≤ ")
    println(h.β)

end



function eliminate_variable(poly_entry, names_variables, var_to_remove)
    poly_exit = eliminate(poly_entry, var_to_remove)
    names_variables_exit = copy(names_variables)
    deleteat!(names_variables_exit, var_to_remove)
    return(poly_exit, names_variables_exit)
end


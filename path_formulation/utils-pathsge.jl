
using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX

includet("../utils/import_utils.jl")







# ========== MASTER PROBLEM
struct MasterProblem
    instance
    model
    gammas
end

struct Gamma
    variable
    path
end



function set_up_master_problem(instance)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    model = Model(CPLEX.Optimizer)

    # ----- Variables
    @variable(model,  0 <= x[v_node in vertices(v_network), s_node in vertices(instance.s_network)] <= 1)


    # ----- Objective
    placement_cost = @expression(model, sum( s_network[s_node][:cost]  * x[v_node, s_node] 
                        for v_node in vertices(v_network), s_node in vertices(s_network) ))
    @objective(model, Min, placement_cost )
    

    # ----- Constraints

    # one substrate node per virtual node
    @constraint(model, [v_node in vertices(v_network)],
        sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1
    )

    # node capacity
    @constraint(model, [s_node in vertices(s_network)],
        sum( x[v_node, s_node] for v_node in vertices(v_network) ) <= s_network[s_node][:cap] 
    )
    
    # one path per v_edge !=> GO FOR >= ???
    @constraint( model, path_selec[v_edge in edges(v_network)],
        0 == 1
    )

    # edge capacity
    @constraint( model, capacity_s_edge[s_edge in edges(s_network)],
        0 <= s_network[src(s_edge), dst(s_edge)][:cap]  
    )

    # start of the path
    @constraint( model, start[v_edge in edges(v_network), s_node in vertices(s_network)],
        0 == x[src(v_edge), s_node]
    )
        
    # terminus
    @constraint( model, terminus[v_edge in edges(v_network), s_node in vertices(instance.s_network)],
        0 == x[dst(v_edge), s_node]
    )


    #= pre-treatment to remove some x - should I use this?
    for v_node in vertices(v_network), s_node in vertices(s_network)
        necessary_bw = degree(v_network, v_node)
        available_bw = sum(s_network_dir[src(s_edge), dst(s_edge)][:cap] for s_edge in get_out_edges(s_node, s_network_dir))
        if necessary_bw > available_bw
            @constraint(model, x[v, s] == 0)
        end
    end
    =#
    add_dumb_columns(instance, model)


    gammas = Dict(v_edge => Vector{Gamma}() for v_edge in edges(v_network))

    return MasterProblem(instance, model, gammas)
end




function add_column(instance, master_problem, v_edge, path)

    model = master_problem.model
    s_network = instance.s_network

    name_col = "γ_$(v_edge)_$(length(master_problem.gammas[v_edge])+1)"
    new_var = @variable(model, base_name=name_col, lower_bound = 0., upper_bound = 1.0)
    push!(master_problem.gammas[v_edge], Gamma(new_var, path))

    set_objective_coefficient(model, new_var, path.cost)
    set_normalized_coefficient(model[:path_selec][v_edge], new_var, 1)
    for s_edge in path.edges
        undirected_s_edge = get_edge(s_network, src(s_edge), dst(s_edge))
        set_normalized_coefficient(model[:capacity_s_edge][undirected_s_edge], new_var, 1)
    end
    set_normalized_coefficient(model[:start][ v_edge, path.src], new_var, 1)
    set_normalized_coefficient(model[:terminus][ v_edge, path.dst], new_var, 1)  
end



function add_dumb_columns(instance, model)

    v_network = instance.v_network
    s_network = instance.s_network
    
    for v_edge in edges(v_network)
        i_dumb = 1
        for s_node in vertices(s_network)
            name_col = "dumb_$(v_edge)_$i_dumb)"
            new_var = @variable(model, base_name=name_col, lower_bound = 0., upper_bound = 1.0)
            set_objective_coefficient(model, new_var, 99999999)
            set_normalized_coefficient(model[:path_selec][v_edge], new_var, 1)
            set_normalized_coefficient(model[:start][ v_edge, s_node], new_var, 1)
            set_normalized_coefficient(model[:terminus][ v_edge, s_node], new_var, 1)  
            i_dumb+=1
        end
    end

end




# ========== DUAL VALUES
struct DualValues
    path_selec::Dict{Edge, Float64}
    capacity_s_edge::Dict{Edge, Float64}
    start::Dict{Edge, Vector{Float64}}
    terminus::Dict{Edge, Vector{Float64}}
end


function get_dual_values(instance, model)

    v_network = instance.v_network
    s_network = instance.s_network

    path_selec= Dict( v_edge=> dual(model[:path_selec][v_edge]) for v_edge in edges(v_network) )
    capacity_s_edge = Dict( s_edge=>dual(model[:capacity_s_edge][s_edge]) for s_edge in edges(s_network) )
    start = Dict( v_edge => [dual( model[:start][v_edge, s_node] ) for s_node in vertices(s_network)] for v_edge in edges(v_network))
    terminus= Dict( v_edge => [dual( model[:terminus][v_edge, s_node] ) for s_node in vertices(s_network)] for v_edge in edges(v_network))

    return DualValues(path_selec, capacity_s_edge, start, terminus)
end


# ========= PRICERS PROBLEMS
function solve_pricer(instance, v_edge, dual_costs)

    # Setting up the mode
    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    model = Model(CPLEX.Optimizer)

    # ----- Variables
    @variable(model,  x_start[s_node in vertices(s_network)], binary=true)
    @variable(model,  x_terminus[s_node in vertices(s_network)], binary=true)
    @variable(model,  y[s_edge in edges(s_network_dir)], binary=true)

    # ----- Objective
    @objective(model, Min, - dual_costs.path_selec[v_edge] +
        sum( (  - dual_costs.start[v_edge][s_node] ) * x_start[s_node] for s_node in vertices(s_network) ) + 
        sum( ( - dual_costs.terminus[v_edge][s_node] ) * x_terminus[s_node] for s_node in vertices(s_network) ) + 
        sum( ( s_network_dir[src(s_edge), dst(s_edge)][:cost] - dual_costs.capacity_s_edge[get_edge(s_network, src(s_edge), dst(s_edge))] ) * y[s_edge] for s_edge in edges(s_network_dir))
    )
    
    
    # ----- Constraints

    # node placement
    @constraint(model, sum( x_start[s_node] for s_node in vertices(s_network)) == 1)
    @constraint(model, sum( x_terminus[s_node] for s_node in vertices(s_network)) == 1)

    # node capacity
    @constraint(model, [s_node in vertices(s_network)],
        x_start[s_node] + x_terminus[s_node] ≤ s_network[s_node][:cap]
    )

    # Flow conservation
    @constraint(model, [s_node in vertices(s_network)],
        sum(  y[s_edge] for s_edge in get_in_edges(s_network_dir, s_node) ) - 
        sum( y[s_edge] for s_edge in get_out_edges(s_network_dir, s_node) ) == 
        x_terminus[s_node] - x_start[s_node] 
    )
    
    # Departure
    @constraint(model, [s_node in vertices(s_network)],
        sum( y[s_edge] for s_edge in get_out_edges(s_network_dir, s_node) ) ≥ x_start[s_node]
    )

    # capacity on nodes
    @constraint(model, [s_node in vertices(s_network)],
        x_start[s_node] + x_terminus[s_node] ≤ 1
    )



    # solve
    set_silent(model)
    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("error! no solution possible...")
        return 
    end

    reduced_cost = objective_value(model)

    if reduced_cost > -0.0001
        return (path_found = nothing, reduced_cost = reduced_cost)
    end

    # Get the solution
    x_start_values = value.(model[:x_start])
    x_terminus_values = value.(model[:x_terminus])
    y_values = value.(model[:y])

    u_start = 0
    u_terminus = 0
    for s_node in  vertices(s_network)
        if x_start_values[s_node] ≥ 0.5
            u_start = s_node
        elseif x_terminus_values[s_node] ≥ 0.5
            u_terminus = s_node
        end
    end

    edges_of_paths = Edge[]
    for  s_edge in edges(s_network_dir)
        if y_values[s_edge] ≥ 0.5
            push!(edges_of_paths, s_edge)
        end
    end

    path = order_path(s_network_dir, edges_of_paths, u_start, u_terminus)

    println("Negative reduced cost! $reduced_cost for $v_edge : $path")


    # return
    return (path_found = path, reduced_cost = reduced_cost)
end





# ===== END HEURISTIC
function make_master_integer(instance, master_problem)

    v_network = instance.v_network
    s_network = instance.s_network
    model = master_problem.model

    for v_node in vertices(v_network), s_node in vertices(s_network)
        set_integer(model[:x][v_node, s_node])
    end



end


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


    model = Model(CPLEX.Optimizer)

    # ----- Setting up the variables constraints etc.
    @variable(model,  0 <= x[vertices(v_network), vertices(instance.s_network)] <= 1)

    # Objective
    placement_cost = @expression(model, sum( s_network[s_node][:cost]  * x[v_node, s_node] 
                for v_node in vertices(v_network) for s_node in vertices(s_network) ))
    @objective(model, Min, placement_cost );


    # constraints

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node ∈ vertices(s_network)) == 1)
    end

    #capacity
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(  x[v_node, s_node] for v_node in vertices(v_network) ) 
                            <= s_network[s_node][:cap] )
    end
    


    #### EDGES

    # one path per v_edge
    @constraint(
        model, 
        path_selec[v_edge in edges(v_network)],
        0 == 1
    );

    # capacity
    @constraint(
        model,
        capacity_s_edge[s_edge in edges(s_network)],
        0 <= s_network[src(s_edge), dst(s_edge)][:cap]  
    );

    # start
    @constraint(
        model, 
        start[v_edge in edges(v_network), s_node in vertices(s_network)],
        0 == x[src(v_edge), s_node]
    );
        
    # terminus
    @constraint(
        model, 
        destination[v_edge in edges(v_network), s_node in vertices(instance.s_network)],
        0 == x[dst(v_edge), s_node]
    );


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

    gammas = Dict()
    for v_edge in edges(v_network)
        gammas[v_edge] = []
    end

    return MasterProblem(instance, model, gammas)
end




function add_column(master_problem, v_edge, path)

    model = master_problem.model
    s_network = instance.s_network

    name_col = "γ_$(src(v_edge))$(dst(v_edge))_$(length(master_problem.gammas[v_edge])+1)"
    new_var = @variable(model, base_name=name_col, lower_bound = 0., upper_bound = 1.0)
    push!(master_problem.gammas[v_edge], Gamma(new_var, path))

    set_objective_coefficient(model, new_var, path.cost)
    set_normalized_coefficient(model[:path_selec][v_edge], new_var, 1)
    for s_edge in path.edges
        set_normalized_coefficient(model[:capacity_s_edge][get_edge(s_network, src(s_edge), dst(s_edge))], new_var, 1)
    end
    set_normalized_coefficient(model[:start][ v_edge, path.src], new_var, 1)
    set_normalized_coefficient(model[:destination][ v_edge, path.dst], new_var, 1)  

end





# ========== DUAL VALUES





# ========= PRICERS PROBLEMS


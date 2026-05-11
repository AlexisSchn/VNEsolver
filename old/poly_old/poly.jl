
using Polyhedra, CDDLib
using XPORTA
includet("../utils/file_reader.jl")
includet("../utils/import_utils.jl")
includet("../resolution/directed/compact/compact_formulation.jl")



function solve_compact_frac_forfun(instance)
    problem = set_up_compact_model_gurobi(instance, true, true, false)
    model = problem.model

    # POUR INSTANCE 2810 1/ : permet de gagner un peu de borne
    #=
    @constraint(model, - model[:x][instance.v_networks[1], 2, 4] - model[:x][instance.v_networks[1], 2, 5] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 6)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 8)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 4)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 6)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 8)]  <= -1)
    @constraint(model, - model[:x][instance.v_networks[1], 1, 4] - model[:x][instance.v_networks[1], 1, 5] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 4)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 6)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 8)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 6)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 8)]  <= -1)
    @constraint(model, - model[:x][instance.v_networks[1], 2, 2] - model[:x][instance.v_networks[1], 2, 3] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 4)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 6)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 8)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 4)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 6)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 8)]  <= -1)
    @constraint(model, - model[:x][instance.v_networks[1], 1, 2] - model[:x][instance.v_networks[1], 1, 3] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 4)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 6)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 8)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 4)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 6)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 8)]  <= -1)
    @constraint(model, - model[:x][instance.v_networks[1], 1, 6] - model[:x][instance.v_networks[1], 1, 7] - model[:x][instance.v_networks[1], 2, 6] - model[:x][instance.v_networks[1], 2, 7] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 4)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 8)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 4)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 8)]  <= -1)
    =#
    #@constraint(model, - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 4)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 6)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 8)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 2)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 4)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 6)] - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 1), get_edge(instance.s_network, 1, 8)]  <= -1)
    
    
    @constraint(model, - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] 
                - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 2, 1)]  
                - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 3, 4)]  
                - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 4, 3)] 
                - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 1, 2)] 
                - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 2, 1)] 
                - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 3, 4)] 
                - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 4, 3)] 
                - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 1, 2)] 
                - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 2, 1)] 
                - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 3, 4)] 
                - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 4, 3)] 

                <= -1)

    # 2=>1, 3=>2, 3=>4
    @constraint(model, - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 2, 1)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 2, 1)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 2, 1)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 3, 2)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 3, 2)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 3, 2)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 3, 4)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 3, 4)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 3, 4)]  
        <= -1)
    # 1=>2, 2=>3, 4=>3
    @constraint(model, - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 1, 2)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 1, 2)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 2, 3)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 2, 3)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 2, 3)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 4, 3)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 4, 3)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 4, 3)]  
        <= -1)

    # 1=>4, 4=>1, 2=>3
    @constraint(model, - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 4)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 1, 4)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 1, 4)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 4, 1)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 4, 1)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 4, 1)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 2, 3)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 2, 3)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 2, 3)]  
        <= -1)

    # 1=>2, 2=>1, 3=>4
    @constraint(model, - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 1, 2)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 1, 2)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 2, 1)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 2, 1)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 2, 1)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 3, 4)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 3, 4)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 3, 4)]  
        <= -1)

    # 3=>2, 4=>1, 1=>4
    @constraint(model, - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 3, 2)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 3, 2)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 3, 2)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 4, 1)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 4, 1)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 4, 1)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 4)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 1, 4)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 1, 4)]  
        <= -1)

    # 2=>3, 4=>1, 1=>4
    @constraint(model, - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 2, 3)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 2, 3)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 2, 3)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 4, 1)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 4, 1)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 4, 1)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 4)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 1, 4)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 1, 4)]  
        <= -1)

    # 4=>3, 3=>4, 1=>2
    @constraint(model, - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 4, 3)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 4, 3)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 4, 3)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 1, 2)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 1, 2)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 1, 2)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 1, 2), get_edge(instance.s_network, 3, 4)] 
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 2, 3), get_edge(instance.s_network, 3, 4)]  
            - model[:y][instance.v_networks[1], get_edge(instance.v_networks[1], 3, 1), get_edge(instance.s_network, 3, 4)]  
        <= -1)

    relax_integrality(model)
    optimize!(model)


    # let's get the solution:

    x_values = value.(problem.x_variables)
    y_values = value.(problem.y_variables)
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            push!(node_placement, [])
            for s_node in vertices(instance.s_network)
                push!(node_placement[v_node], x_values[v_network, v_node, s_node])
            end
        end
        edge_routing = Dict()
        total_val = 0
        for v_edge in edges(v_network)
            edge_routing[v_edge] = Dict()
            val = 0
            for s_edge in edges(instance.s_network)
                val += y_values[v_network, v_edge, s_edge]
                edge_routing[v_edge][s_edge] = y_values[v_network, v_edge, s_edge]
            end
            total_val += val
            println("For edge $v_edge we have $val")
        end
        println("Total val : $total_val")
        m = MappingCompactFractional(v_network, instance.s_network, node_placement, edge_routing)
        print(m)
    end


end


function get_vrep_undir(instance)

    instance_dir = get_directed_instance(instance)
    problem = set_up_compact_model_gurobi(instance_dir, true, true, true)
    model = problem.model
    set_optimizer_attribute(model, "PoolSearchMode", 2)
    set_optimizer_attribute(model, "PoolSolutions", 1000)


    optimize!(model)
    #solution_summary(model)
    println("There are $(result_count(model)) solution in the model")
    sols = Vector{Vector{Int64}}()
    for i in 1:result_count(model)
        cur_sol = Vector{Int64}()
        for v_node in vertices(instance.v_networks[1])
            for s_node in vertices(instance_dir.s_network)
                push!(cur_sol, round.(value.(model[:x][instance_dir.v_networks[1], v_node, s_node]; result = i)))     
            end
        end

        for v_edge in edges(instance.v_networks[1])
            for s_edge in edges(instance_dir.s_network)
                push!(cur_sol, round.(value.(model[:y][instance_dir.v_networks[1], v_edge, s_edge]; result = i)))     
            end
        end

        push!(sols, cur_sol)
    end

    names = []
    for v_node in vertices(instance.v_networks[1])
        for s_node in vertices(instance.s_network)
            push!(names, "x_" * string(v_node) * "_" * string(s_node))
        end
    end
    
    for v_edge in edges(instance.v_networks[1])
        for s_edge in edges(instance_dir.s_network)
            push!(names, "y_" * string(src(v_edge)) * string(dst(v_edge)) * "_" * string(src(s_edge)) * string(dst(s_edge)))
        end
    end

    return sols, names
end


function get_vrep_dir(instance)

    problem = set_up_compact_model_gurobi(instance, true, true, false)
    model = problem.model
    set_optimizer_attribute(model, "PoolSearchMode", 2)
    set_optimizer_attribute(model, "PoolSolutions", 10000)


    optimize!(model)
    #solution_summary(model)
    println("There are $(result_count(model)) solution in the model")
    sols = Vector{Vector{Int64}}()
    for i in 1:result_count(model)
        cur_sol = Vector{Int64}()
        for v_node in vertices(instance.v_networks[1])
            for s_node in vertices(instance.s_network)
                push!(cur_sol, round.(value.(model[:x][instance.v_networks[1], v_node, s_node]; result = i)))     
            end
        end

        for v_edge in edges(instance.v_networks[1])
            for s_edge in edges(instance.s_network)
                push!(cur_sol, round.(value.(model[:y][instance.v_networks[1], v_edge, s_edge]; result = i)))     
            end
        end

        push!(sols, cur_sol)
    end

    names = []
    for v_node in vertices(instance.v_networks[1])
        for s_node in vertices(instance.s_network)
            push!(names, "x_" * string(v_node) * "_" * string(s_node))
        end
    end
    
    for v_edge in edges(instance.v_networks[1])
        for s_edge in edges(instance.s_network)
            push!(names, "y_" * string(src(v_edge)) * string(dst(v_edge)) * "_" * string(src(s_edge)) * string(dst(s_edge)))
        end
    end

    return sols, names
end



# takes an instance (with one vr !), returns all the solutions (points) and the names of variables
function get_v_rep_convexhull(instance, time_limit=100, sol_max=1000)

    
    # builds the model and find all optimal points (a bit obscure in cplex, i found this online)
    compact_model = set_up_compact_model(instance)
    model = compact_model.model
    #set_time_limit_sec(model, time_limit)

    set_silent(model)
    optimize!(model)
    #set_optimizer_attribute(model, "CPX_PARAM_DISPLAY", 0)

    set_optimizer_attribute(model, "CPX_PARAM_SOLNPOOLAGAP", 5) # we accept +5 compared to the best
    set_optimizer_attribute(model, "CPX_PARAM_SOLNPOOLINTENSITY", 2) # Not super aggressive populate, maybe could try 3 ?
    set_optimizer_attribute(model, "CPX_PARAM_POPULATELIM", 10)
    
    backend_model = unsafe_backend(model);
    env = backend_model.env;
    lp = backend_model.lp;
    
    CPLEX.CPXpopulate(env, lp);
    
    N_results = CPLEX.CPXgetsolnpoolnumsolns(env, lp)
    if N_results == sol_max
        println("The number of solution max has been reached. Some solutions might not have been compute.")
    end

    # get correct names
    v_network = instance.v_networks[1]
    names_variables_dic = Dict()
    
    for v_node in vertices(v_network)
        for s_node in vertices(instance.s_network)
            names_variables_dic[CPLEX.column(backend_model, compact_model.x_variables[v_network, v_node, s_node].index)] = "x_" * string(v_node) * "_" * string(s_node)
        end
    end
    
    for v_edge in edges(v_network)
        for s_edge in edges(instance.s_network)
            names_variables_dic[CPLEX.column(backend_model, compact_model.y_variables[v_network, v_edge, s_edge].index)] = "y_" * string(src(v_edge)) * string(dst(v_edge)) * "_" * string(src(s_edge)) * string(dst(s_edge))
        end
    end
    
    names_variables = []
    for i_var in 1:length(names_variables_dic)
        push!(names_variables, names_variables_dic[i_var])
    end

    # get all solutions vectors
    all_sols = Vector{Vector{Int64}}()
    for sol in 0:N_results-1
        current_sol = Vector{Int64}()
        for i_variable in 0:length(names_variables)-1
            x = Ref{Cdouble}()  ## Reference to the `Term` variable value
            CPLEX.CPXgetsolnpoolx(env, lp, sol, x, i_variable, i_variable)
            push!(current_sol, convert.(Int64, round.(x[])))
        end
        push!(all_sols, current_sol)
    end

    println("\n\nThere are " * string(N_results) * " solutions")


    
    return(all_sols, names_variables)
end



# takes an instance (with one vr !), returns all the solutions (points) and the names of variables
function get_v_rep_convexhull_undir(instance_undir, time_limit=100, sol_max=1000)

    instance = get_directed_instance(instance_undir)

    # builds the model and find all optimal points (a bit obscure in cplex, i found this online)
    compact_model = set_up_compact_model(instance)
    model = compact_model.model
    #set_time_limit_sec(model, time_limit)

    set_silent(model)
    optimize!(model)
    #set_optimizer_attribute(model, "CPX_PARAM_DISPLAY", 0)

    set_optimizer_attribute(model, "CPX_PARAM_SOLNPOOLAGAP", 100000000.0)
    set_optimizer_attribute(model, "CPX_PARAM_SOLNPOOLINTENSITY", 4)
    set_optimizer_attribute(model, "CPX_PARAM_POPULATELIM", sol_max)
    
    backend_model = unsafe_backend(model);
    env = backend_model.env;
    lp = backend_model.lp;
    
    CPLEX.CPXpopulate(env, lp);
    
    N_results = CPLEX.CPXgetsolnpoolnumsolns(env, lp)
    if N_results == sol_max
        println("The number of solution max has been reached. Some solutions might not have been compute.")
    end

    # get correct names
    v_network = instance.v_networks[1]
    
    names_variables = []
    sols = []
    for i in 1:N_results
        push!(sols, [])
    end

    for v_node in vertices(v_network)
        for s_node in vertices(instance.s_network)
            name = "x_" * string(v_node) * "_" * string(s_node)
            push!(names_variables, name)

            col = CPLEX.column(backend_model, compact_model.x_variables[v_network, v_node, s_node].index)
            println("col $col : " * name)
            for (i_sol, sol) in enumerate(sols)
                x = Ref{Cdouble}() 
                CPLEX.CPXgetsolnpoolx(env, lp, i_sol, x, col, col)
                if (x[]) > 0.5
                    println("in sol $i_sol")
                end
                if (x[]) < - 0.5
                    println(" WHAT THE FUUUUUUUUUCK")
                end
                push!(sol, convert.(Int64, round.(x)))

            end
        end
    end
    
    for v_edge in edges(instance_undir.v_networks[1])
        for s_edge in edges(instance_undir.s_network)
            #check da stuff
        end
    end
    

    println("\n\nThere are " * string(N_results) * " solutions")

    for (i_sol, sol) in enumerate(all_sols)
        println("Printing sol $i_sol...")
        for (i_var, val) in enumerate(sol)            
            if val > 0.5
                println("$(names_variables[i_var])")
            end
            
        end
    end

    return(sols, names_variables)
end





function print_polytope(poly, names_variables)
    hr = hrep(poly);
    println("There are " * string(length(names_variables)) * " variables.\n")
    println("There are " * string(length(hyperplanes(hr))) * " hyperplanes")
    for h in hyperplanes(hr)
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                print("+ ")
                print(floor(Int, h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                print("\t")
            elseif (h.a[i_var] < -0.0001)
                print("- ")
                print(floor(Int, -h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                print("\t")
            else
                print("\t\t")
            end
        end
        print(" = ")
        println(floor(Int,h.β))
    end

    println("\n\n There are " * string(length(halfspaces(hr))) * " halfspaces")

    for h in halfspaces(hr)
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                print("+ ")
                print(floor(Int, h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                print("\t")
            elseif (h.a[i_var] < -0.0001)
                print("- ")
                print(floor(Int, -h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                print("\t")
            else
                print("\t\t")
            end
        end
        print(" ≤ ")
        println(floor(Int,h.β))
    end

end


function print_polytope_simpler(poly, names_variables, print_index=false)
    hr = hrep(poly);
    println("There are " * string(length(names_variables)) * " variables.\n")
    println("There are " * string(length(hyperplanes(hr))) * " hyperplanes")
    for h in hyperplanes(hr)
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                print("+ ")
                print(floor(Int, h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                if print_index
                    print("[" * string(i_var) * "]")
                end
                print("\t")
            elseif (h.a[i_var] < -0.0001)
                print("- ")
                print(floor(Int, -h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                if print_index
                    print("[" * string(i_var) * "]")
                end
                print("\t")
            end
        end
        print(" = ")
        println(floor(Int,h.β))
    end

    println("\n\n There are " * string(length(halfspaces(hr))) * " halfspaces")

    for h in halfspaces(hr)
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                print("+ ")
                print(floor(Int, h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                if print_index
                    print("[" * string(i_var) * "]")
                end
                print("\t")
            elseif (h.a[i_var] < -0.0001)
                print("- ")
                print(floor(Int, -h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                if print_index
                    print("[" * string(i_var) * "]")
                end
                print("\t")
            end
        end
        print(" ≤ ")
        println(floor(Int,h.β))
    end



end

function print_polytope_simpler_smaller(poly, names_variables, print_index=false)
    hr = hrep(poly);
    println("There are " * string(length(names_variables)) * " variables.\n")
    println("There are " * string(length(hyperplanes(hr))) * " hyperplanes")
    for h in hyperplanes(hr)
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                print("+ ")
                print(floor(Int, h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                if print_index
                    print("[" * string(i_var) * "]")
                end
                print("\t")
            elseif (h.a[i_var] < -0.0001)
                print("- ")
                print(floor(Int, -h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                if print_index
                    print("[" * string(i_var) * "]")
                end
                print("\t")
            end
        end
        print(" = ")
        println(floor(Int,h.β))
    end

    println("\n\n There are " * string(length(halfspaces(hr))) * " halfspaces")
    #=
    for h in halfspaces(hr)
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                print("+ ")
                print(floor(Int, h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                if print_index
                    print("[" * string(i_var) * "]")
                end
                print("\t")
            elseif (h.a[i_var] < -0.0001)
                print("- ")
                print(floor(Int, -h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                if print_index
                    print("[" * string(i_var) * "]")
                end
                print("\t")
            end
        end
        print(" ≤ ")
        println(floor(Int,h.β))
    end

    =#

end


function print_variable(names_variables)
    for i_var in 1:length(names_variables)
        println(string(i_var) * " : " * names_variables[i_var])
    end
end


function eliminate_variables(poly_entry, names_variables, to_elim)
    poly_exit = eliminate(poly_entry, to_elim)
    names_variables_exit = copy(names_variables)
    for i_truc in 1:length(to_elim)
        deleteat!(names_variables_exit, to_elim[length(to_elim) - i_truc + 1])
    end
    return(poly_exit, names_variables_exit)
end


function do_everything(instance)
    all_sols, names_variables = get_vrep_dir(instance);
    v_rep = vrep(all_sols);
    poly = polyhedron(v_rep, CDDLib.Library(:float));
    return poly, names_variables
end

function elimate_and_do_everything(poly_entry, names_variables, to_elim)
    poly_project = eliminate(poly_entry, to_elim)
    names_variables_project = copy(names_variables)
    for i_truc in 1:length(to_elim)
        deleteat!(names_variables_project, to_elim[length(to_elim) - i_truc + 1])
    end
    return poly_project, names_variables_project
end

function do_everything_and_print(instance)
    all_sols, names_variables = get_vrep_dir(instance);
    v_rep = vrep(all_sols);
    poly = polyhedron(v_rep, CDDLib.Library(:float));
    h_rep = hrep(poly);
    print_polytope(h_rep, names_variables);
end



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


function print_constraints(poly, instance)

    hr = hrep(poly);
    #=
    for h in hyperplanes(hr)
        print("@constraint(model, ")
        
        i_var = 1
        for (i_vn, v_network) in enumerate(instance.v_networks)
            for v_node in vertices(v_network)
                for s_node in vertices(instance.s_network)
                    if (h.a[i_var] > 0.01)
                        print("+ model[:x][instance.v_networks[$i_vn], $v_node, $s_node] ")
                    elseif (h.a[i_var] < -0.01)
                        print("- model[:x][instance.v_networks[$i_vn], $v_node, $s_node] ")
                    end
                    i_var += 1
                end
            end

            for v_edge in edges(v_network)
                for s_edge in edges(instance.s_network)
                    if (h.a[i_var] > 0.01)
                        print("+ model[:y][instance.v_networks[$i_vn], get_edge(instance.v_networks[$i_vn], $(src(v_edge)), $(dst(v_edge))), get_edge(instance.s_network, $(src(s_edge)), $(dst(s_edge)))] ")
                    elseif (h.a[i_var] < -0.01)
                        print("- model[:y][instance.v_networks[$i_vn], get_edge(instance.v_networks[$i_vn], $(src(v_edge)), $(dst(v_edge))), get_edge(instance.s_network, $(src(s_edge)), $(dst(s_edge)))] ")
                    end
                    i_var += 1
                end
            end
            print(" == ")
            print(floor(Int,h.β))
            println(")")
        end

    end
    =#

    for h in halfspaces(hr)
        if (floor(Int,h.β)< -0.5)
            print("@constraint(model, ")
            
            i_var = 1
            for (i_vn, v_network) in enumerate(instance.v_networks)
                for v_node in vertices(v_network)
                    for s_node in vertices(instance.s_network)
                        if (h.a[i_var] > 0.01)
                            print("+ model[:x][instance.v_networks[$i_vn], $v_node, $s_node] ")
                        elseif (h.a[i_var] < -0.01)
                            print("- model[:x][instance.v_networks[$i_vn], $v_node, $s_node] ")
                        end
                        i_var += 1
                    end
                end

                for v_edge in edges(v_network)
                    for s_edge in edges(instance.s_network)
                        if (h.a[i_var] > 0.01)
                            print("+ model[:y][instance.v_networks[$i_vn], get_edge(instance.v_networks[$i_vn], $(src(v_edge)), $(dst(v_edge))), get_edge(instance.s_network, $(src(s_edge)), $(dst(s_edge)))] ")
                        elseif (h.a[i_var] < -0.01)
                            print("- model[:y][instance.v_networks[$i_vn], get_edge(instance.v_networks[$i_vn], $(src(v_edge)), $(dst(v_edge))), get_edge(instance.s_network, $(src(s_edge)), $(dst(s_edge)))] ")
                        end
                        i_var += 1
                    end
                end
                print(" <= ")
                print(floor(Int,h.β))
                println(")")
            end
        end

    end

end


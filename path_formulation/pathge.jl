

includet("../utils/import_utils.jl")
includet("utils-pathsge.jl")

using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX



function solve_pathge(instance)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir
    master_problem = set_up_master_problem(instance)


    # add original paths: those that are just one edge long...
    for s_edge in edges(instance.s_network_dir)
        cost = s_network_dir[src(s_edge), dst(s_edge)][:cost]
        if s_network[src(s_edge)][:cap] > 0 && s_network[dst(s_edge)][:cap] > 0
            s_path = Path(src(s_edge), dst(s_edge), [s_edge], cost)
            for v_edge in edges(v_network)
                add_column(instance, master_problem, v_edge, s_path)
            end
        end
    end
    set_silent(master_problem.model)

    model = master_problem.model
    optimize!(model)
    
    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Unfeasible!")
        return (10e9)
    end
    
    println("Well the objective value is : $(objective_value(model))")

    dual_costs = get_dual_values(instance, model)

    keep_on = true
    while keep_on
        keep_on = false
        # pricers
        for v_edge in edges(v_network)
            result = solve_pricer(instance, v_edge, dual_costs)
            if result[:reduced_cost] ≤ -0.001
                add_column(instance, master_problem, v_edge, result[:path_found])
                keep_on = true
            end
        end
    
        # rmp solving
        optimize!(model)
        println("CG value: $(objective_value(model))")
        dual_costs = get_dual_values(instance, model)
    end

    println("\n\n\nFinal value: $(objective_value(model))")


    #println("Now, let's do it integer:")
    #unset_silent(master_problem.model)
    #make_master_integer(instance, master_problem)
    #optimize!(master_problem.model)

    return round(objective_value(model); digits=3)
end

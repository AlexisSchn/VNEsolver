

includet("../utils/import_utils.jl")
includet("utils-pathsge.jl")

using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX



function solve_pathge(instance)

    v_network = instance.v_network
    s_network = instance.s_network

    master_problem = set_up_master_problem(instance)


    # add original paths: those that are just one edge long...
    for s_edge in edges(instance.s_network_dir)
        cost = instance.s_network_dir[src(s_edge), dst(s_edge)][:cost]
        s_path = Path(src(s_edge), dst(s_edge), [s_edge], cost)
        for v_edge in edges(v_network)
            add_column(master_problem, v_edge, s_path)
        end
    end

    optimize!(master_problem.model)

    println("Well the objective value is : $(objective_value(master_problem.model))")

    #=
    for s_edge in edges(s_network)
        cap_used = 0
        for v_edge in edges(v_network)
            for gamma in master_problem.gammas[v_edge]
                for s_edge_path in gamma.path.edges
                    if get_edge(s_network, src(s_edge_path), dst(s_edge_path)) == s_edge
                        cap_used += value.(gamma.variable)
                    end
                end
            end
        end
        println("For $s_edge, cap used: $cap_used, cap: $(s_network[src(s_edge), dst(s_edge)][:cap])")
    end
    =#
    # solve the master problem for that shit.


    #=
    pricers = Dict()
    for v_edge in edges(v_network)
        pricers[v_edge] = set_up_pricer(v_network)
    end


    while keep_on

        keep_on = false

        optimize!(master_problem.model)

        print("CG value: $(objective_value(model))")
        for v_edge in edges(g_network)
            column, value = update_and_solve_pricer(pricers[v_edge], dual_values)
            if value < -0.001
                add_column(master_problem, column)
                keep_on = true
            end

        end
    
    end
    =#

    return

end

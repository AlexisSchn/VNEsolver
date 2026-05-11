using Revise, JuMP, CPLEX, Gurobi, Metis
using Graphs, MetaGraphsNext

includet("../../../utils/import_utils.jl")
includet("../../../utils/visu.jl")



function test_cycle(instance)

    println("Hello there !")
    
    cycles_basis_v_network = cycle_basis(instance.v_network)

    v_node = 13

    # get triangle in which v_node is in
    triangle_with_vnode = []
    for cycle in cycles_basis_v_network
        if v_node ∈ cycle && length(cycle) == 3
            println("New v triangle : $cycle")
            push!(triangle_with_vnode, cycle)
        end
    end

    # check if s_node is in a triangle. Otherwise, remove it.

    cycles_basis_s_network = cycle_basis(instance.s_network)

    v_node = 9

    # get triangle in which v_node is in
    triangle_with_snode = []
    for cycle in cycles_basis_s_network
        if v_node ∈ cycle && length(cycle) == 3
            println("New s triangle : $cycle")
            push!(triangle_with_snode, cycle)
        end
    end


    visu_graph(instance.s_network)
 

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
    
    
    
    # Outgoing edges cap: pretty stupid but useful
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
    #println("We get this to delete: $i")
    
end




function add_some_cycle_ineqs(instance, model)


    s_network = instance.s_network
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir

    # take some 5 v_nodes
    cycles_basis_v_network = cycle_basis(v_network)
    triangles_basis_v_network = []
    for cycle in cycles_basis_v_network
        if length(cycle) == 3
            push!(triangles_basis_v_network, cycle)
        end
    end

    # get the substrate nodes that are not in any triangle
    cycles_basis_s_network = cycle_basis(s_network)
    triangle_basis_s_network = []
    for cycle in cycles_basis_s_network
        if length(cycle) == 3
            push!(triangle_basis_s_network, cycle)
        end
    end

    s_nodes_in_no_triangles = []
    for s_node in vertices(s_network)
        no_triangle = true
        for triangle in triangle_basis_s_network
            if s_node ∈ triangle
                no_triangle = false
            end
        end
        if no_triangle
            push!(s_nodes_in_no_triangles, s_node)
        end
    end

    println("Here are the s_node in no triangles : $s_nodes_in_no_triangles")

    nb_constraints = 0
    # Let's gongue
    for triangle in triangles_basis_v_network
        v_edges_triangle = [] #get em..

        for v_node_1 in triangle
            for v_node_2 in triangle
                if has_edge(v_network, v_node_1, v_node_2)
                    v_edge = get_edge(v_network, v_node_1, v_node_2)
                    if v_edge ∉ v_edges_triangle
                        push!(v_edges_triangle, v_edge)
                    end
                end
            end
        end
        println("Edges of triangle $triangle: $v_edges_triangle")

        
        for v_node in triangle
            @constraint(model, sum(model[:y][v_edge, s_edge] for v_edge in v_edges_triangle for s_edge in edges(s_network_dir)) 
                >=
                3 + sum( model[:x][v_node, s_node] for s_node in s_nodes_in_no_triangles)    
            )
            nb_constraints +=1
        end
        
    end

    println("We have added $nb_constraints constraints ")

end




function solve_compact_cycles(instance, time_solver = 5, stay_silent=false)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir


    model = Model(CPLEX.Optimizer)
    set_up_problem(instance, model)
    add_some_cycle_ineqs(instance, model)

    #set_optimizer_attribute(model, "CPXPARAM_Emphasis_MIP", 1)

    #relax_integrality(problem.model)
    set_time_limit_sec(model, time_solver)
    if stay_silent
        set_silent(model)
    else
        print("Starting solving... ")
    end
    
    optimize!(model)

    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("Infeasible or unfinished: $status")
        return
    end


end



function solve_compact_cycles_and_stars(instance, time_solver = 5, stay_silent=false)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir


    model = Model(CPLEX.Optimizer)
    set_up_problem(instance, model)
    add_some_cycle_ineqs(instance, model)
    add_star_constraints(instance, model)
    #set_optimizer_attribute(model, "CPXPARAM_Emphasis_MIP", 1)

    #relax_integrality(problem.model)
    set_time_limit_sec(model, time_solver)
    if stay_silent
        set_silent(model)
    else
        print("Starting solving... ")
    end
    
    optimize!(model)

    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("Infeasible or unfinished: $status")
        return
    end


end



function solve_compact(instance, time_solver = 5, stay_silent=false)
    

    model = Model(CPLEX.Optimizer)
    set_up_problem(instance, model)

    set_time_limit_sec(model, time_solver)
    if stay_silent
        set_silent(model)
    else
        print("Starting solving... ")
    end
    
    optimize!(model)

    #=
    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("Infeasible or unfinished: $status")
        return
    end

    x_values = value.(model[:x])
    y_values = value.(model[:y])

    if !stay_silent
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
    end
    =#

end


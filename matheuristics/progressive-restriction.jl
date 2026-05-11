
using Revise, JuMP, CPLEX, Gurobi
includet("../../../utils/import_utils.jl")


struct Compact_Formulation
    model
    x_variables
    y_variables
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




    #=

    #----------  Local constraints 1
    diff_deg=1
    δ2_edges = Dict()
    for s_node in vertices(s_network_dir)
        edges_in_distance = δ2(s_network_dir, s_node)
        #println("Edge in δ2 of $s_node : $edges_in_distance")
        δ2_edges[s_node] = edges_in_distance
    end

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
                if v_deg - s_deg > diff_deg
                    cons = @constraint(model, 
                        trucmuche >= (v_deg - s_deg) * model[:x][v_node, s_node]
                    )
                end
            end
        end
    end



    #----------  Local constraints 2
    diff_deg=1
    for (v_node, v_deg) in dense_v_nodes
        v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
        for s_node in vertices(s_network_dir)
            s_deg = degree(instance.s_network, s_node) 
            if v_deg - s_deg > diff_deg

                # Here we need to be extra careful about the pseudo orientation of the virtual network
                trucmuche = @expression(model, 0.)
                for v_edge in v_edges_incident
                    if src(v_edge) == v_node
                        trucmuche += sum(model[:y][v_edge, s_edge] for neighbor in neighbors(s_network_dir, s_node) for s_edge in get_out_edges(s_network_dir, neighbor) )
                    else
                        trucmuche += sum(model[:y][v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  for neighbor in neighbors(s_network_dir, s_node) for s_edge in get_out_edges(s_network_dir, neighbor) )
                    end
                end
                
                @constraint(model, 
                    trucmuche >= (v_deg - s_deg) * model[:x][v_node, s_node] + sum(v_deg * model[:x][v_node, s_neighbor] for s_neighbor in neighbors(s_network_dir, s_node))
                )
            end
        end
    end
    =#



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


function set_up_problem_placement_restrict(instance, model, placement_restriction)

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
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1) ### A RAJOUTER ABSOLUMNET ! (ou alors tu ne créés pas les autres variables).
        # sinon il est pas au courant qu'il peut "supprimer" les x pas valables...

        if v_node in keys(placement_restriction)
            @constraint(model, sum(x[v_node, s_node] for s_node in placement_restriction[v_node]) == 1)
        end
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

    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end

    #--------------- infeasibility constraints: if there is just not enough capacity on edges yknow...
    nb_var_less = 0
    for v_node in vertices(v_network)
        for s_node in vertices(s_network)
            v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
            necessary_bw = sum(v_network[src(v_edge), dst(v_edge)][:dem] for v_edge in v_edges_incident)

            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident)
            if necessary_bw > available_bw
                nb_var_less += 1
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end
    
    
end




function set_up_problem_placement_restrict_modified(instance, model, placement_restriction)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    ### Variables
    @variable(model, x[vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);

    

    ### Objective

    ####" !!! We are going to modify the cost, why not, I can't help to push v_node elsewhere for now so that's all I can try :////
    presence_s_node = []
    for s_node in vertices(s_network)
        pres = 0
        for v_node in keys(placement_restriction)
            if s_node ∈ placement_restriction[v_node]
                pres+=1
            end
        end
        push!(presence_s_node, pres)
    end


    println("Presence s nodes : $(presence_s_node)")

    #println("Modified cost :")

    modified_cost = Dict()
    for s_node in vertices(s_network)
        modified_cost[s_node] = Dict()
        for v_node in vertices(s_network)
            cost = 0
            change = true
            if v_node ∈ keys(placement_restriction)
                if s_node ∈ placement_restriction[v_node]
                    cost = s_network[s_node][:cost]
                    change = false
                end
            end
            if change
                cost = s_network[s_node][:cost] * (1 + presence_s_node[s_node]*3)
            end
            modified_cost[s_node][v_node] = cost

        end
        #println("For $s_node : $(modified_cost[s_node])")
    end


    placement_cost = @expression(model, sum( modified_cost[s_node][v_node] * v_network[v_node][:dem] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);




    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1) ### A RAJOUTER ABSOLUMNET ! (ou alors tu ne créés pas les autres variables).
        # sinon il est pas au courant qu'il peut "supprimer" les x pas valables...

        if v_node in keys(placement_restriction)
            @constraint(model, sum(x[v_node, s_node] for s_node in placement_restriction[v_node]) == 1)
        end
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

    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end

    #--------------- infeasibility constraints: if there is just not enough capacity on edges yknow...
    nb_var_less = 0
    for v_node in vertices(v_network)
        for s_node in vertices(s_network)
            v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
            necessary_bw = sum(v_network[src(v_edge), dst(v_edge)][:dem] for v_edge in v_edges_incident)

            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident)
            if necessary_bw > available_bw
                nb_var_less += 1
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end
    
    
end


function solve_heuristic_relax(instance)

    v_network = instance.v_network
    s_network = instance.s_network

    nb_s_nodes= 15

    nodes_to_restrict_next = []


    nb_neighbors = []
    nb_neighbors_done = []
    for v_node in vertices(v_network)
        push!(nb_neighbors, length(neighbors(v_network, v_node)))
        push!(nb_neighbors_done, 0)
    end


    restrictions = Dict()
    restricted_nodes = []

    # we simply takes the most dense node and make it start here
    sorted_degree = sort(collect(vertices(v_network)), by=v_node -> -degree(v_network, v_node))
    nodes_to_restrict_next = sorted_degree[1:1]


    iter=0
    while length(restricted_nodes) != length(vertices(v_network))
        iter+=1
        current_model = Model(CPLEX.Optimizer)
        println("\n---------------\nIter $iter\n---------------")
        #println("Current placement restriction: $(restrictions)")
        set_up_problem_placement_restrict_modified(instance, current_model, restrictions) 
        add_star_constraints(instance, current_model)
        relax_integrality(current_model)
        set_silent(current_model)
    
        optimize!(current_model)

        if !has_values(current_model)
            println("Infeasible model !")
            return
        end
        #println("Solving time: $(solve_time(current_model)) ")
        println("Value LP: $(objective_value(current_model)), solved in $(solve_time(current_model))")
        
        x_values = value.(current_model[:x])


        for v_node in nodes_to_restrict_next
            selec = []
            println("For virtual node $(v_node)")
            for i in 1:nb_s_nodes
                #=
                vals = [(x_values[v_node, s_node] * 100
                        + sum(x_values[v_node, neighbor] * 5 for neighbor in neighbors(s_network, s_node); init=0.0)
                        - sum(x_values[v_node, neighbor] *10 for neighbor in intersect(neighbors(s_network, s_node), selec); init=0.0)
                        - sum(x_values[already_chosen_vnode, s_node] * 40 /length(restricted_nodes) for already_chosen_vnode in restricted_nodes; init=0.0)) 
                            for s_node in vertices(s_network)]
                =#

                vals = [(x_values[v_node, s_node] * 100
                        + sum(x_values[v_node, neighbor] * 5 for neighbor in neighbors(s_network, s_node); init=0.0)) 
                            for s_node in vertices(s_network)]
                sorted = sortperm(vals, rev=true)
                i_truc = 1
                stop = false 
                while !stop
                    if sorted[i_truc] ∉ selec
                        push!(selec, sorted[i_truc])
                        stop = true
                    end
                    i_truc += 1
                end
            end
            restrictions[v_node] = selec
            println("Placement restriction for node $v_node: $selec")

            push!(restricted_nodes, v_node)

        end


                # -------- NEXT NODES SELECTION: quite tricky ?
        # idea : 1) most central left, 2) most neighbors 3) less frac 4) all neighbors ?
        next_nodes_choice = 3
        if length(restricted_nodes) != length(vertices(v_network))

            if next_nodes_choice == 1
                ratio = Dict()
                for v_node in vertices(v_network)
                    if v_node ∉ restricted_nodes
                        #ratio[v_node] = -nb_neighbors_done[v_node] / nb_neighbors[v_node] # with less neighbor already fixed first
                        #ratio[v_node] = nb_neighbors_done[v_node] / nb_neighbors[v_node] # with most neighbor already fixed first
                        ratio[v_node] = degree(v_network, v_node) # with highest degree node
                    end
                end
                sorted_keys = sort(collect(keys(ratio)), by=x -> -ratio[x])
                nodes_to_restrict_next = [sorted_keys[1]]
            elseif next_nodes_choice == 2

                nodes_to_restrict_next= []
                for v_node in restricted_nodes
                    for neighbor in neighbors(v_network, v_node)
                        if neighbor ∉ restricted_nodes
                            if neighbor ∉ nodes_to_restrict_next
                                push!(nodes_to_restrict_next, neighbor)
                            end
                        end
                    end
                end

            elseif next_nodes_choice == 3 # On en prend qu'un, le plus dense parmis les voisins... ?
                possible_nodes = []
                for v_node in restricted_nodes
                    for neighbor in neighbors(v_network, v_node)
                        if neighbor ∉ restricted_nodes
                            if neighbor ∉ possible_nodes
                                push!(possible_nodes, neighbor)
                            end
                        end
                    end
                end
                ratio = Dict()
                for v_node in possible_nodes
                    #ratio[v_node] = -nb_neighbors_done[v_node] / nb_neighbors[v_node] # with less neighbor already fixed first
                    #ratio[v_node] = nb_neighbors_done[v_node] / nb_neighbors[v_node] # with most neighbor already fixed first
                    ratio[v_node] = degree(v_network, v_node) # with highest degree node
                end
                sorted_keys = sort(collect(keys(ratio)), by=x -> -ratio[x])
                nodes_to_restrict_next = [sorted_keys[1]]

            end

            println("Next nodes to restrict: $nodes_to_restrict_next")


        end

    end

    # model improved
    model_restrict = Model(CPLEX.Optimizer)
    set_up_problem_placement_restrict(instance, model_restrict, restrictions)
    print("Starting solving model_restrict ")
    optimize!(model_restrict)

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
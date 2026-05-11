
using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf

#general
includet("../utils/import_utils.jl")

# utils colge
includet("utils/master_problem.jl")
includet("utils/graph_decomposition_overlapping.jl")
includet("utils/column_generation.jl")
includet("../utils/partition-graph.jl")

# init
includet("init/greedy_easy.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")




function node_decompo(instance)

    println("Starting...")
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    # ======= SETTING UP THE DECOMPOSITION ======= #

    # AUTOMATIC PARTITION


    v_node_partitionning = []
    for v_node in vertices(v_network)
        push!(v_node_partitionning, [v_node])
    end


    println("Node partitionning: $v_node_partitionning")





    vn_decompo = set_up_decompo_overlapping(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    
    # === COLUMN GENERATION === #

    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")



    # column generation!
    return column_generation(instance, vn_decompo, master_problem)


end



function edge_decompo(instance)

    println("Starting...")
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    # ======= SETTING UP THE DECOMPOSITION ======= #

    # AUTOMATIC PARTITION


    v_node_partitionning = []
    for v_edge in edges(v_network)
        push!(v_node_partitionning, [src(v_edge), dst(v_edge)])
    end


    println("Node partitionning: $v_node_partitionning")





    vn_decompo = set_up_decompo_overlapping(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    
    # === COLUMN GENERATION === #

    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")



    # generating first columns. For edge...
    #=
    for s_edge in edges(s_network_dir)

        for v_subgraph in vn_subgraphs
            node_placement = [src(s_edge), dst(s_edge)]
            cost_mapping = s_network_dir[src(s_edge)][:cost] + s_network_dir[dst(s_edge)][:cost] + s_network_dir[src(s_edge), dst(s_edge)][:cost]
            edge_routing = Dict( collect(edges(v_subgraph.graph))[1] => Path(src(s_edge), dst(s_edge), [s_edge], s_network_dir[src(s_edge), dst(s_edge)][:cost]))
            sub_mapping = Mapping(v_subgraph.graph, s_network_dir, node_placement, edge_routing)
            add_column(master_problem, instance, vn_decompo, v_subgraph, sub_mapping, cost_mapping)
        end

    end
    =#

    # column generation!
    return column_generation(instance, vn_decompo, master_problem)


    # ======= END HEURISTIC STUFF ======= #

    #basic_heuristic(instance, vn_decompo, master_problem, 900)

    return 
end




function stars_partition(instance)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    subgraphs = []
    copy_v_network = copy(v_network.graph)
    keep_on = true

    real_indices = collect(1:nv(v_network)) 
    # Graphs.jl is completly stupid when it comes to removing a node. It swaps it with the last node and then removes it. be careful... 

    nodes_in_no_subgraphs = collect(1:nv(v_network))
    v_node_partitionning = []
    while keep_on

        # Get node with max degree
        node_with_max_degree = argmax(degree(copy_v_network))

        keep_on = true
        if degree(copy_v_network, node_with_max_degree) == 0
            break
        end

        #println("most central node: $node_with_max_degree, aka $(real_indices[node_with_max_degree])")

        new_part = [ real_indices[node_with_max_degree] ]

        edges_subg = []
        for neigh in neighbors(copy_v_network, node_with_max_degree)
            #println("Neighbor: $neigh, aka $(real_indices[neigh])")
            push!(new_part, real_indices[neigh])
            push!(edges_subg, get_edge(v_network, real_indices[node_with_max_degree], real_indices[neigh]))
        end

        push!(v_node_partitionning, new_part)
        
        #push!(subgraphs, Dict("nodes"=>new_part, "edges"=>edges_subg))

        for v_node in new_part
            if v_node ∈ nodes_in_no_subgraphs
                filter!(x -> x != v_node, nodes_in_no_subgraphs)
            end
        end

        rem_vertex!(copy_v_network, node_with_max_degree)
        real_indices[node_with_max_degree] = real_indices[length(real_indices)]
        deleteat!(real_indices, length(real_indices))

        
        #println("Real indices: $real_indices")

    end


    #println("Some nodes left? $nodes_in_no_subgraphs")
    #println("Node partitionning: $v_node_partitionning")

    for v_node_left in nodes_in_no_subgraphs
        push!(v_node_partitionning, [v_node_left])
    end

    println("Node partitionning: $v_node_partitionning")





    vn_decompo = set_up_decompo_overlapping(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges: $(vn_decompo.v_edges_master)")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    

    # === COLUMN GENERATION === #

    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")

    # column generation!
    column_generation(instance, vn_decompo, master_problem)

end




# A bit less stars, with a bit less overlapping nodes?
# The center of a star can not be in another star.
function stars_partition_2(instance)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    copy_v_network = copy(v_network.graph)
    keep_on = true

    real_indices = collect(1:nv(v_network)) 
    # Graphs.jl is completly stupid when it comes to removing a node. It swaps it with the last node and then removes it. be careful... 

    nodes_in_no_subgraphs = collect(1:nv(v_network))
    v_node_partitionning = []
    possible_centers = collect(1:nv(v_network))

    while keep_on

        # Get node with max degree
        nodes_degrees = [degree(copy_v_network, v_node) for v_node in vertices(copy_v_network)]
        node_sorted_degree = sortperm(nodes_degrees, rev=true)
        center = 0
        for node in node_sorted_degree
            if real_indices[node] ∈ possible_centers
                center = node
                break
            end
        end
        #println("The center will be $center")
        if center == 0
            #println("No center found?")
            break
        end
        if degree(copy_v_network, center) == 0
            #println("Degree too low: time to stop.")
            break
        end
        #println("New center! $center")
        #println("most central node: $node_with_max_degree, aka $(real_indices[node_with_max_degree])")

        new_part = [ real_indices[center] ]

        edges_subg = []
        for neigh in neighbors(copy_v_network, center)
            #println("Neighbor: $neigh, aka $(real_indices[neigh])")
            push!(new_part, real_indices[neigh])
            push!(edges_subg, get_edge(v_network, real_indices[center], real_indices[neigh]))
        end

        # Nodes from this star can not be centers anymore.
        for v_node in new_part
            filter!(x -> x != v_node, possible_centers)
        end

        push!(v_node_partitionning, new_part)
        
        for v_node in new_part
            if v_node ∈ nodes_in_no_subgraphs
                filter!(x -> x != v_node, nodes_in_no_subgraphs)
            end
        end

        rem_vertex!(copy_v_network, center)
        real_indices[center] = real_indices[length(real_indices)]
        deleteat!(real_indices, length(real_indices))

        
        #println("Nodes that can be center: $possible_centers")
        #println("Nodes that are still in the graph: $real_indices")
    end


    #println("Some nodes left? $nodes_in_no_subgraphs")
    #println("Node partitionning: $v_node_partitionning")

    for v_node_left in nodes_in_no_subgraphs
        push!(v_node_partitionning, [v_node_left])
    end

    println("Node partitionning: $v_node_partitionning")





    vn_decompo = set_up_decompo_overlapping(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges: $(vn_decompo.v_edges_master)")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    

    # === COLUMN GENERATION === #

    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")

    # column generation!
    column_generation(instance, vn_decompo, master_problem)


end


# Here, looking at real stars:
# Subgraphs are real stars, I don't take edges that are "wheel edges".
# More edges in master problem
function stars_overlapping_decompo(instance)
    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    subgraphs = []
    copy_v_network = copy(v_network.graph)
    keep_on = true

    real_indices = collect(1:nv(v_network)) 
    # Graphs.jl is completly stupid when it comes to removing a node. It swaps it with the last node and then removes it. be careful... 

    nodes_in_no_subgraphs = collect(1:nv(v_network))
    v_node_partitionning = []
    possible_centers = collect(1:nv(v_network))

    while keep_on

        # Get node with max degree
        nodes_degrees = [degree(copy_v_network, v_node) for v_node in vertices(copy_v_network)]
        node_sorted_degree = sortperm(nodes_degrees, rev=true)
        center = 0
        for node in node_sorted_degree
            if real_indices[node] ∈ possible_centers
                center = node
                break
            end
        end
        #println("The center will be $center")
        if center == 0
            #println("No center found?")
            break
        end
        if degree(copy_v_network, center) == 0
            #println("Degree too low: time to stop.")
            break
        end

        #println("most central node: $node_with_max_degree, aka $(real_indices[node_with_max_degree])")

        new_part = [ real_indices[center] ]

        edges_subg = []
        for neigh in neighbors(copy_v_network, center)
            #println("Neighbor: $neigh, aka $(real_indices[neigh])")
            push!(new_part, real_indices[neigh])
            push!(edges_subg, get_edge(v_network, real_indices[center], real_indices[neigh]))
        end

        # Nodes from this star can not be centers anymore.
        for v_node in new_part
            filter!(x -> x != v_node, possible_centers)
        end

        
        push!(subgraphs, Dict("nodes"=>new_part, "edges"=>edges_subg))

        for v_node in new_part
            if v_node ∈ nodes_in_no_subgraphs
                filter!(x -> x != v_node, nodes_in_no_subgraphs)
            end
        end

        rem_vertex!(copy_v_network, center)
        real_indices[center] = real_indices[length(real_indices)]
        deleteat!(real_indices, length(real_indices))

        
        #println("Real indices: $real_indices")

    end


    println("Some nodes left? $nodes_in_no_subgraphs")
    #println("Node partitionning: $v_node_partitionning")




    vn_decompo = set_up_decompo_overlapping_more_info(instance, subgraphs)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges: $(vn_decompo.v_edges_master)")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    

    # === COLUMN GENERATION === #
    
    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")

    # column generation!
    return column_generation(instance, vn_decompo, master_problem)
    #column_generation_greedy(instance, vn_decompo, master_problem)

end




# No overlappin nodes!
function stars_strict_decompo(instance)
    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    subgraphs = []
    copy_v_network = copy(v_network.graph)
    keep_on = true

    real_indices = collect(1:nv(v_network)) 
    # Graphs.jl is completly stupid when it comes to removing a node. It swaps it with the last node and then removes it. be careful... 

    nodes_in_no_subgraphs = collect(1:nv(v_network))
    v_node_partitionning = []
    possible_centers = collect(1:nv(v_network))

    while keep_on

        # Get node with max degree
        nodes_degrees = [degree(copy_v_network, v_node) for v_node in vertices(copy_v_network)]
        node_sorted_degree = sortperm(nodes_degrees, rev=true)
        center = 0
        for node in node_sorted_degree
            if real_indices[node] ∈ possible_centers
                center = node
                break
            end
        end
        #println("The center will be $center")
        if center == 0
            #println("No center found?")
            break
        end
        if degree(copy_v_network, center) == 0
            #println("Degree too low: time to stop.")
            break
        end

        #println("most central node: $node_with_max_degree, aka $(real_indices[node_with_max_degree])")

        new_part = [ real_indices[center] ]

        edges_subg = []
        for neigh in neighbors(copy_v_network, center)
            #println("Neighbor: $neigh, aka $(real_indices[neigh])")
            push!(new_part, real_indices[neigh])
            push!(edges_subg, get_edge(v_network, real_indices[center], real_indices[neigh]))
        end

        # Nodes from this star can not be centers anymore.
        for v_node in new_part
            filter!(x -> x != v_node, possible_centers)
        end

        
        push!(subgraphs, Dict("nodes"=>new_part, "edges"=>edges_subg))

        for v_node in new_part
            if v_node ∈ nodes_in_no_subgraphs
                filter!(x -> x != v_node, nodes_in_no_subgraphs)
            end
        end


        # remove first the leafs of the star, because of how it's done in Graphs.jl
        while length(neighbors(copy_v_network, center)) > 0
            #println("Well now, the neighboring of $center has $(length(neighbors(copy_v_network, center))) nodes..")
            neigh = neighbors(copy_v_network, center)[1]
            rem_vertex!(copy_v_network, neigh)

            real_indices[neigh] = real_indices[length(real_indices)]
            #println("Real indices: $(length(real_indices))")
            if center == length(real_indices)
                center = neigh
                #println("Well it's time to change, now the center is $center")
            end
            deleteat!(real_indices, length(real_indices))


        end

        rem_vertex!(copy_v_network, center)
        real_indices[center] = real_indices[length(real_indices)]
        deleteat!(real_indices, length(real_indices))

    end


    println("Some nodes left? $nodes_in_no_subgraphs")
    #println("Node partitionning: $v_node_partitionning")

    # Add them to the decompo?
    for v_node in nodes_in_no_subgraphs
        push!(subgraphs, Dict("nodes"=>[v_node], "edges"=>[]))
    end



    vn_decompo = set_up_decompo_overlapping_more_info(instance, subgraphs)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges: $(vn_decompo.v_edges_master)")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    

    # === COLUMN GENERATION === #
    
    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")

    # column generation!
    #column_generation_greedy(instance, vn_decompo, master_problem)
    return column_generation(instance, vn_decompo, master_problem)

end



# paths !
# overlapping paths!
# maybe triangle and cycles ? and the reminder are paths ?
# It would be nice to use cycle base.


function path_overlapping_decompo(instance)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    subgraphs = []
    copy_v_network = copy(v_network.graph)
    keep_on = true


    # While the the shortest path in the residual network is longer than 1, you keep on
    while keep_on
        
        shortest_paths = floyd_warshall_shortest_paths(copy_v_network)
        #println("Yo dists: $(shortest_paths.dists)")
        for i in 1:length(shortest_paths.dists)
            if  shortest_paths.dists[i] > 50.
                shortest_paths.dists[i] = 0.
            end
        end
        #println("After? dists: $(shortest_paths.dists)")

        if maximum(shortest_paths.dists) < 0.5
            break
        end
        
        couple = argmax(shortest_paths.dists)
        print("Wow the max length in the residual v network is $couple with $(maximum(shortest_paths.dists))")



        src = couple[1]
        dst = couple[2]
        nodes_of_path = [dst]
        edges_of_path = []
        v = dst
        while v != src
            u = shortest_paths.parents[src, v]
            push!(nodes_of_path, u)
            push!(edges_of_path, get_edge(v_network, u, v))

            rem_edge!(copy_v_network, get_edge(copy_v_network, u, v))

            v = shortest_paths.parents[src, v]
        end

        println(" The path is: $nodes_of_path")

        push!(subgraphs, Dict("nodes"=>nodes_of_path, "edges"=>edges_of_path))

    end


    vn_decompo = set_up_decompo_overlapping_more_info(instance, subgraphs)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges: $(vn_decompo.v_edges_master)")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    

    # === COLUMN GENERATION === #
    
    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")

    # column generation!
    return column_generation(instance, vn_decompo, master_problem)



end



function path_strict_decompo(instance)


    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    subgraphs = []
    copy_v_network = copy(v_network.graph)
    keep_on = true
    nodes_in_no_subgraphs = collect(1:nv(v_network))

    # While the the shortest path in the residual network is longer than 1, you keep on
    while keep_on
        
        shortest_paths = floyd_warshall_shortest_paths(copy_v_network)
        #println("Yo dists: $(shortest_paths.dists)")
        for i in 1:length(shortest_paths.dists)
            if  shortest_paths.dists[i] > 50.
                shortest_paths.dists[i] = 0.
            end
        end
        #println("After? dists: $(shortest_paths.dists)")

        if maximum(shortest_paths.dists) < 0.5
            break
        end
        
        couple = argmax(shortest_paths.dists)
        print("Wow the max length in the residual v network is $couple with $(maximum(shortest_paths.dists))")



        src = couple[1]
        dst = couple[2]
        nodes_of_path = [dst]
        edges_of_path = []
        v = dst
        while v != src
            u = shortest_paths.parents[src, v]
            push!(nodes_of_path, u)
            push!(edges_of_path, get_edge(v_network, u, v))
            v = shortest_paths.parents[src, v]
        end

        println(" The path is: $nodes_of_path")

        push!(subgraphs, Dict("nodes"=>nodes_of_path, "edges"=>edges_of_path))

        for v_node in nodes_of_path
            neighs = copy(neighbors(copy_v_network, v_node))
            for neigh in neighs
                rem_edge!(copy_v_network, get_edge(copy_v_network, v_node, neigh))
            end
        end

        nodes_in_no_subgraphs = setdiff(nodes_in_no_subgraphs, nodes_of_path)

    end

    println("Some nodes left? $nodes_in_no_subgraphs")
    for v_node in nodes_in_no_subgraphs
        push!(subgraphs, Dict("nodes"=>[v_node], "edges"=>[]))
    end

    vn_decompo = set_up_decompo_overlapping_more_info(instance, subgraphs)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges: $(vn_decompo.v_edges_master)")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    

    # === COLUMN GENERATION === #
    
    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")

    # column generation!
    return column_generation(instance, vn_decompo, master_problem)



end



function cycle_overlapping_decompo(instance)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    subgraphs = []
    copy_v_network = copy(v_network.graph)
    keep_on = true
    nodes_in_no_subgraphs = collect(1:nv(v_network))


    while keep_on
        
        if length(cycle_basis(copy_v_network)) == 0
            break
        end
        nb_cycle = argmin([length(bas) for bas in cycle_basis(copy_v_network)])
        cycle = cycle_basis(copy_v_network)[nb_cycle]
        println("Best cycle: $cycle")

        edges_of_cycle= []

        for i_node in 1:length(cycle)-1
            push!(edges_of_cycle, get_edge(v_network, cycle[i_node], cycle[i_node+1]))
        end
        push!(edges_of_cycle, get_edge(v_network, cycle[length(cycle)], cycle[1]))

        push!(subgraphs, Dict("nodes"=>cycle, "edges"=>edges_of_cycle))

        for v_edge in edges_of_cycle
            rem_edge!(copy_v_network, v_edge)
        end

        nodes_in_no_subgraphs = setdiff(nodes_in_no_subgraphs, cycle)

    end


    println("Some nodes left? $nodes_in_no_subgraphs")
    for v_node in nodes_in_no_subgraphs
        push!(subgraphs, Dict("nodes"=>[v_node], "edges"=>[]))
    end

    vn_decompo = set_up_decompo_overlapping_more_info(instance, subgraphs)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges: $(vn_decompo.v_edges_master)")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    

    # === COLUMN GENERATION === #
    
    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")

    # column generation!
    return column_generation(instance, vn_decompo, master_problem)



end



function cycle_strict_decompo(instance)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    subgraphs = []
    copy_v_network = copy(v_network.graph)
    keep_on = true
    nodes_in_no_subgraphs = collect(1:nv(v_network))


    while keep_on
        
        if length(cycle_basis(copy_v_network)) == 0
            break
        end
        nb_cycle = argmin([length(bas) for bas in cycle_basis(copy_v_network)])
        cycle = cycle_basis(copy_v_network)[nb_cycle]
        println("Best cycle: $cycle")

        edges_of_cycle= []

        for i_node in 1:length(cycle)-1
            push!(edges_of_cycle, get_edge(v_network, cycle[i_node], cycle[i_node+1]))
        end
        push!(edges_of_cycle, get_edge(v_network, cycle[length(cycle)], cycle[1]))

        push!(subgraphs, Dict("nodes"=>cycle, "edges"=>edges_of_cycle))


        for v_node in cycle
            neighs = copy(neighbors(copy_v_network, v_node))
            for neigh in neighs
                rem_edge!(copy_v_network, get_edge(copy_v_network, v_node, neigh))
            end
        end

        nodes_in_no_subgraphs = setdiff(nodes_in_no_subgraphs, cycle)

    end


    println("Some nodes left? $nodes_in_no_subgraphs")
    for v_node in nodes_in_no_subgraphs
        push!(subgraphs, Dict("nodes"=>[v_node], "edges"=>[]))
    end

    vn_decompo = set_up_decompo_overlapping_more_info(instance, subgraphs)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges: $(vn_decompo.v_edges_master)")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    

    # === COLUMN GENERATION === #
    
    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")

    # column generation!
    return column_generation(instance, vn_decompo, master_problem)



end







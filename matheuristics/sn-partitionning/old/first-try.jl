using Revise, JuMP, CPLEX, Gurobi, Metis
using Graphs, MetaGraphsNext

includet("../../../utils/import_utils.jl")
includet("../../../utils/visu.jl")

includet("../../compact_undir.jl")


function solve_partitionning(instance)

    v_network = instance.v_network
    s_network = instance.s_network

    println("Let's take the s_network and partition it")

    nb_clusters = Int(floor(nv(s_network) /  nv(v_subgraph.graph)))
    println("We need $(nb_s_subgraph) substrate subgraphs for $v_subgraph.")
    partition = Metis.partition(s_network.graph, nb_clusters)

    part = Dict(i=>Vector{Int64}() for i in 1:nb_clusters)
    for i in 1:nv(s_network)
        push!(part[partition[i]], i)
    end

    for a in keys(part)
        println("There are $(length(part[a])) nodes in $a : $(part[a])")
    end

    println("What do you think ? Visualize")
    visu_partitioning(s_network.graph, partition)

    
    for i_cluster in 1:nb_clusters
        println("Now let's cook for cluster $i_cluster:\n")
        sub_s_network = my_induced_subgraph(s_network, part[i_cluster], "sn_$i_cluster")
        subinstance = Instance_Undir_VNE_1s(v_network, sub_s_network)
        solve_compact_1vn(subinstance, 30)
        println("\n\n\n")
    end
    
end





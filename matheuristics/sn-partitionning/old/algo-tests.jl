using Revise, JuMP, CPLEX, Gurobi, Metis
using Graphs, MetaGraphsNext

includet("../../../utils/import_utils.jl")
includet("../../../utils/visu.jl")



function metis_test(g, nb_clusters=3, type_algo=2)

    println("Let's take the s_network and partition it for $nb_clusters clusters")
    if type_algo ==1
        partition = Metis.partition(g, nb_clusters; alg=:KWAY)
    else
        partition = Metis.partition(g, nb_clusters; alg=:RECURSIVE)
    end

    part = Dict(i=>Vector{Int64}() for i in 1:nb_clusters)
    for i in 1:nv(g)
        push!(part[partition[i]], i)
    end

    for a in keys(part)
        println("There are $(length(part[a])) nodes in $a : $(part[a])")
    end

    println("What do you think ? Visualize")
    visu_partitioning(g, partition)

    
end



function clique_perco_test(instance, perco_number = 0)

    s_network = instance.s_network
    partition = clique_percolation(s_network.graph, k=2)
    println(partition)

end

function label_propagation_test(instance)

    s_network = instance.s_network
    partition = label_propagation(s_network, 1000)[1]

    visu_partitioning(s_network.graph, partition)
end
using Base.Libc

# Path to the KaHIP shared library
const KAHIP_LIB = "/home/alexiss/KaHIP/build/libkahip.so"  # Update with the actual path


using Revise, JuMP, CPLEX, Gurobi, Metis
using Graphs, MetaGraphsNext

includet("../../../utils/import_utils.jl")
includet("../../../utils/visu.jl")




function partition_kahip(g, nb_clusters, imbalance)
    xadj, adjncy = graph_to_csr(g)

    for i in 1:length(adjncy)
        adjncy[i] = adjncy[i] - 1
    end

    # adapt it for c, because they start indexing at 0 lol

    edgecut, part = simplified_kaffpa(xadj, adjncy, nb_clusters, imbalance)

    # adapt it for julia, because we start indexing at 1 lol
    for i in 1:length(part)
        part[i] = part[i] + 1
    end

    return part
end



function do_everything_kahip(g, nb_clusters)
    xadj, adjncy = graph_to_csr(g)
    println("Xadj: $xadj")
    println("Adjency: $adjncy")


    for i in 1:length(adjncy)
        adjncy[i] = adjncy[i] - 1
    end

    # adapt it for c, because they start indexing at 0 lol

    imbalance = 0.20
    edgecut, part = simplified_kaffpa(xadj, adjncy, nb_clusters, imbalance)

    # adapt it for julia, because we start indexing at 1 lol
    for i in 1:length(part)
        part[i] = part[i] + 1
    end

    println("Edge cut: $edgecut")
    println("Partition: $part")

    println("Let's vizualize. What do you think ?")
    visu_partitioning(g, part)
end



function graph_to_csr(g::Graph)
    n = nv(g)  # Number of vertices

    # Build xadj and adjncy arrays
    xadj = Int[]
    adjncy = Int[]  # This will hold all adjacency lists concatenated

    push!(xadj, 0)
    for i in 1:n
        append!(adjncy, neighbors(g, i))
        push!(xadj, length(adjncy))
    end

    return xadj, adjncy
end



function simplified_kaffpa(xadj, adjncy, nparts, imbalance)
    n = Int32(length(xadj)-1)  # Number of vertices

    suppress_output = true  # Suppress output
    seed = Int32(42)  # Random seed
    mode = Int32(0)  # Mode (e.g., FAST)
    
    nparts = Int32(nparts)
    imbalance = Float64(imbalance)

    xadj = Vector{Int32}(xadj)
    adjncy = Vector{Int32}(adjncy)

    edgecut, part = kaffpa(n, xadj, adjncy, nparts, imbalance, suppress_output, seed, mode)

    #println("Edge cut: $edgecut")
    #println("Partition: $part")

    return(edgecut, part)

end


function kaffpa(n::Int32,
                xadj::Vector{Int32},
                adjncy::Vector{Int32}, 
                nparts::Int32, 
                imbalance::Float64, 
                suppress_output::Bool, 
                seed::Int32, 
                mode::Int32)
    
    # Prepare output parameters
    edgecut = Ref{Int32}(0)  # Edge cut result
    part = Vector{Int32}(undef, n)  # Partition array to be filled by kaffpa

    # Call kaffpa
    ccall(
        (:kaffpa, KAHIP_LIB), 
        Cvoid,  # Return type
        (   
            Ptr{Cint},      #n
            Ptr{Cint},      #vwgt
            Ptr{Cint},      #xadj
            Ptr{Cint},      #edge weights
            Ptr{Cint},      #Adjacency
            Ptr{Cint},      # n parts
            Ptr{Cdouble},   #inbalance
            Bool, 
            Cint, 
            Cint, 
            Ptr{Cint}, 
            Ptr{Cint}
        ),
        Ref(n),
        C_NULL, 
        xadj, 
        C_NULL, 
        adjncy, 
        Ref(nparts), 
        Ref(imbalance), 
        suppress_output, 
        seed, 
        mode, 
        edgecut, 
        part
    )

    # Return results
    return edgecut[], part
end

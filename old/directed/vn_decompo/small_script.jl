using Revise
includet("../../utils/file_reader.jl")
includet("../../utils/some_data.jl")
includet("../../utils/utils_graph.jl")
includet("../../resolution/compact/compact_formulation.jl")
includet("../../resolution/path/path_generation.jl")
includet("../../utils/visu.jl")
includet("vn_decompo.jl")

node_partitionning = [[[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16, 17]]]

dir = "../../instances/tests/vndemo/"
for file in readdir(dir)
    println("\n\n\n\n\nCurrent instance: " * dir * file)
    instance = get_instance_from_folder(dir*file)
    println("SN with: " * string(length(vertices(instance.s_network))) * " nodes and " * string(length(edges(instance.s_network))) * " edges.")

    println("\n Starting PLNE solving with 10s:")
    solve_compact_integer(instance, true, true, 10, true);

    time1 = time()
    println("\n Starting PLNE solving with 100s:")
    solve_compact_integer(instance, true, true, 100, true);
    println("in "  * string(time1 - time()) * "seconds.")

    println("\n\nStarting CG:")
    vn_decompo(instance, node_partitionning);
end
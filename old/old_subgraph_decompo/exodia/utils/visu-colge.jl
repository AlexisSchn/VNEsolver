using Graphs, MetaGraphsNext
using GraphRecipes
using Plots
using Colors
using NetworkLayout


function layoutaaahhh(g)



end

function visu_graph(g)
    w = []
    for i_node in 1:nv(g)
        if i_node < 10
            push!(w, 10)
        else
            push!(w, 1)
        end
    end


    
    coord = spring(g, iterations=500)
    p = graphplot(g,
        x=[coord[i][1]*10 for i in 1:nv(g)],
        y=[coord[i][2]*10 for i in 1:nv(g)],
        node_weights = w,
        names=string.(1:nv(g)),
        curvature_scalar=0.01, 
        node_size = 8
    )
    display(p) 

end


function visu_repartition(g, repartition)
    w = []
    for i_node in 1:nv(g)
        if i_node < 10
            push!(w, 100)
        elseif i_node < 100
            push!(w, 10)
        else
            push!(w, 1)
        end
    end



    colors = distinguishable_colors(maximum(repartition), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    #colors = distinguishable_colors(maximum(repartition))

    marker_cols = []
    for i_node in 1:nv(g)
        push!(marker_cols, colors[repartition[i_node]])
    end

    coord = spring(g, iterations=500)
    p = graphplot(g,
        x=[coord[i][1]*0.20 for i in 1:nv(g)],
        y=[coord[i][2]*0.20 for i in 1:nv(g)],
        node_weights = w,
        names=string.(1:nv(g)),
        curvature_scalar=0.01,
        markercolor=marker_cols,
        node_size = 0.5,
        dpi=300
    )
    display(p) 

    return p
end

function visu_partitionning(g, partitionning)

    repartition = []

    for node in vertices(g)
        for (i_part, part) in enumerate(partitionning)
            if node ∈ part
                push!(repartition, i_part)
            end
        end     
    end

    p = visu_repartition(g, repartition)
    return p
end


function visu_added_nodes(g, nodes, added)
    w = []
    for i_node in 1:nv(g)
        if i_node < 10
            push!(w, 100)
        elseif i_node < 100
            push!(w, 10)
        else
            push!(w, 1)
        end
    end

    colors = distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)

    marker_cols = []
    for i_node in 1:nv(g)
        if i_node ∈ added
            push!(marker_cols, colors[2])
        elseif i_node ∈ nodes
            push!(marker_cols, colors[3])
        else
            push!(marker_cols, colors[1])
        end
    end

    coord = spring(g, iterations=500)
    p = graphplot(g,
        x=[coord[i][1]*0.20 for i in 1:nv(g)],
        y=[coord[i][2]*0.20 for i in 1:nv(g)],
        node_weights=w,
        names=string.(1:nv(g)),
        markercolor=marker_cols,
        curvature_scalar=0.01,
        node_size = 0.4,
        dpi=300 
    )

    return p
end
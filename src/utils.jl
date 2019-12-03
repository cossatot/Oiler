#module Utils

#using ..BlockRotations

include("./block_rotations.jl")

function diagonalize_matrices(matrices)

    rowz = [size(m, 1) for m in matrices]
    colz = [size(m, 2) for m in matrices]

    n_rows = sum(rowz)
    n_cols = sum(colz)

    big_mat = zeros(n_rows, n_cols)

    i_row = 1
    i_col = 1

    for (i, mat) in enumerate(matrices)
        row_end = i_row + rowz[i] - 1
        col_end = i_col + colz[i] - 1

        big_mat[i_row:row_end, i_col:col_end] = mat

        i_row = row_end + 1
        i_col = col_end + 1
    end
    return big_mat
end


function make_digraph_from_vels(vels::VelocityVectorSphere...)
#function make_digraph_from_vels(vels)

    vel_graph = Dict{String, Array{String}}()

    for vel in vels
        if haskey(vel_graph, vel.fix)
            if vel.mov in vel_graph[vel.fix]
            # pass
            else
                push!(vel_graph[vel.fix], vel.mov)
            end
        else
            vel_graph[vel.fix] = [vel.mov]
        end
    end
    vel_graph
end



#function make_digraph_from_poles(poles::Union{EulerPoleCart, EulerPoleSphere})
#
#    pole_graph = SimpleDiGraph()
#
#    for pole in poles
#        if pole.fix != 
#    end

#end
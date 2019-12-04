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
    vel_array = collect(vels)
    make_digraph_from_vels(vel_array)
end


"""
    make_digraph_from_vels(vels)

Makes a directed graph from a group of velocity vectors. Each velocity vector
(of type VelocityVectorSphere) has 'fix' and 'mov' fields that describe the
velocity of a point or block (mov) relative to another point or block (fix). The
directed graph maps out these relationships to determine constraints on block
motions, in particular ensuring closure of velocity circuits. Edges are not
duplicated, so if multiple velocities constrain relative motion between two
blocks, only one edge is given in the resulting graph.

The graph is in an adjacency list format.

# Arguments

# Returns


"""
function make_digraph_from_vels(vels::Array{VelocityVectorSphere})
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


function make_ugraph_from_digraph(digraph::Dict)
    ug = Dict{String, Array{String}}()

    for (fix, movs) in digraph
        for mov in movs
            if !haskey(ug, fix)
                ug[fix] = [mov]
            elseif !(mov in ug[fix])
                push!(ug[fix], mov)
            end
            if !haskey(ug, mov)
                ug[mov] = [fix]
            elseif !(fix in ug[mov])
                push!(ug[mov], fix)
            end
        end
    end
    ug
end


#function make_digraph_from_poles(poles::Union{EulerPoleCart, EulerPoleSphere})
#
#    pole_graph = SimpleDiGraph()
#
#    for pole in poles
#        if pole.fix != 
#    end

#end
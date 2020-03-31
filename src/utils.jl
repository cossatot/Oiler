module Utils

export predict_vels_from_poles, find_vel_cycles, diagonalize_matrices, 
    random_sample_vel_groups, build_Vc_from_vel_samples, get_gnss_vels,
    get_coords_from_vel_array

using ..Oiler: VelocityVectorSphere, PoleCart, PoleSphere, build_PvGb_from_vels,
        build_vel_column_from_vels, add_poles, pole_sphere_to_cart


using DataFrames
using SparseArrays
using LinearAlgebra


function diagonalize_matrices(matrices)

    rowz = [size(m, 1) for m in matrices]
    colz = [size(m, 2) for m in matrices]

    n_rows = sum(rowz)
    n_cols = sum(colz)

    big_mat = spzeros(n_rows, n_cols)

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
    vel_graph = Dict{String,Array{String}}()

    for vel in vels
        if haskey(vel_graph, vel.fix)
            if !(vel.mov in vel_graph[vel.fix])
                push!(vel_graph[vel.fix], vel.mov)
            end
        else
            vel_graph[vel.fix] = [vel.mov]
        end
    end
    vel_graph
end


function make_digraph_from_poles(poles::Union{Array{PoleCart},Array{PoleSphere}})
    pole_graph = Dict{String,Array{String}}()

    for pole in poles
        if haskey(pole_graph, pole.fix)
            if !(pole.mov in pole_graph[pole.fix])
                push!(pole_graph[pole.fix], pole.mov)
            end
        else
            pole_graph[pole.fix] = [pole.mov]
        end
    end
    pole_graph
end


function make_digraph_from_tuples(tups)
    graph = Dict()

    for (t1, t2) in tups
        if haskey(graph, t1)
            if !(t2 in graph[t1])
                push!(graph[t1], t2)
            end
        else
            graph[t1] = [t2]
        end
    end
    graph
end

function make_ugraph_from_digraph(digraph::Dict)
    ug = Dict{String,Array{String}}()

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


function find_tricycles(graph::Dict; reduce::Bool = true)

    all_tris = []

    for (from, tos) in graph
        for to in tos
            if to != from
                for next in graph[to]
                    if next != from
                        for fin in graph[next]
                            if fin == from
                                path = [from, to, next]
                                push!(all_tris, path)
                            end
                        end
                    end
                end
            end
        end
    end

    if reduce == true
        sets = []
        unique_tris = []

        for tri in all_tris
            if !(Set(tri) in sets)
                push!(sets, Set(tri))
                push!(unique_tris, tri)
            end
        end
        return unique_tris
    else 
        return all_tris
    end
end


function get_cycle_inds(vel_group_keys, cycle_tup)

    v_ind = findall(x->x == cycle_tup, vel_group_keys)
    if v_ind != []
        res = Dict("ind" => v_ind[1], "val" => 1.)
    else
        cycle_tup = Base.reverse(cycle_tup)
        v_ind = findall(x->x == (cycle_tup), vel_group_keys)
        if v_ind != []
            res = Dict("ind" => v_ind[1], "val" => -1.)
        else
            println("v apparently not in cycle")
        end
    end
    res
end


function find_vel_cycles(vel_group_keys)

    cycle_inds = Dict()
    dg = make_digraph_from_tuples(vel_group_keys)
    tricycles = find_tricycles(make_ugraph_from_digraph(dg))

    for (i, tri) in enumerate(tricycles)
        cycle_inds[i] = Dict()
        for ctup in [(tri[1], tri[2]),(tri[2], tri[3]),(tri[3], tri[1])]
            cycle_inds[i][ctup] = get_cycle_inds(vel_group_keys, ctup)
        end
    end
    cycle_inds
end


"""
    flat(x)
Flattens nested arrays of arrays
"""
flat(x, y = vcat(x...)) = x == y ? x : flat(y)


"""
    find_shortest_path(graph, start, stop)

Finds the shortest path between vertices `start` and `stop` given the
edges in the (directed, unweighted) `graph` using a breadth-first-search
strategy.

# Arguments
- `graph`: A dictionary with keys for each vertex, and values of arrays of
  all the vertices linked to the key vertex.

- `start`: Vertex to start from.

- `stop`: Vertex to stop on.

# Returns
Array of vertices representing the shortest path in order from `start` to 
`stop`.
"""
function find_shortest_path(graph::Dict{String,Array{String}}, 
    start::String, stop::String)

    dist = Dict{String,Any}(start => [start])
    q = [start]
    while length(q) > 0
        at = popfirst!(q)
        for next in graph[at]
            if !haskey(dist, next)
                dist[next] = [dist[at], next]
                push!(q, next)
            end
        end
    end
    flat(dist[stop])
end


function get_pole_path(poles::Array{PoleCart}, path::Array{String})

    steps = length(path) - 1
    pole_path = Array{PoleCart,1}(undef, steps)


    for i in 1:steps
        place = path[i]
        next = path[i + 1]

        for pole in poles
            if pole.fix == place && pole.mov == next
                pole_path[i] = pole
            elseif pole.fix == next && pole.mov == place
                pole_path[i] = -pole
            end
        end
    end
    pole_path
end


function get_path_euler_pole(poles::Array{PoleCart,1}, fix::String,
mov::String)
    
    if fix == mov
        final_pole = PoleCart(x = 0., y = 0., z = 0., fix = fix, mov = mov)
    else
        vel_dg = make_digraph_from_poles(poles)
        vel_ug = make_ugraph_from_digraph(vel_dg)

        path = find_shortest_path(vel_ug, fix, mov)

        pole_path = get_pole_path(poles, path)

        final_pole = add_poles(pole_path)
    end
    final_pole
end


function
predict_vels_from_poles(block_things::Dict{String,AbstractArray},
    poles::Array{PoleCart,1})

    pole_list = [get_path_euler_pole(poles, fix, mov) for (fix, mov) in
    block_things["keys"]]
    
    pole_vec = reduce(vcat, [[p.x; p.y; p.z] for p in pole_list])

    V_pred = block_things["PvGb"] * pole_vec

    Vn_pred = V_pred[1:3:end]
    Ve_pred = V_pred[2:3:end]
    # Vu_pred = V_pred[3:3:end]

    (Ve_pred, Vn_pred)
end

function predict_vels_from_poles(block_things::Dict{String,AbstractArray},
    poles::Array{PoleSphere})

    cpoles = [pole_sphere_to_cart(pole) for pole in poles]

    predict_vels_from_poles(block_things, cpoles)
end


"""
    get_gnss_vels(vel_groups)

Collects all of the GNSS-derived velocities from the velocity groups, with
some ancillary metadata.

# Arguments

# Returns

"""
function get_gnss_vels(vel_groups)
    gnss_vels = []
    row_set_num = 0
    for (i, group) in enumerate(values(vel_groups))
        col_idx = 3 * (i - 1) + 1
        for vel in group
            row_set_num += 1
            if vel.vel_type == "GNSS"
                row_idx = 3 * (row_set_num - 1) + 1
                vel_idx = [row_idx:row_idx + 2, col_idx:col_idx + 2]
                vd = Dict()
                vd["vel"] = vel
                vd["idx"] = vel_idx
                push!(gnss_vels, vd)
            end # if
        end # for
    end # for
    gnss_vels
end


"""
    get_coords_from_vel_array(vels)

Returns (lons, lats) as arrays of floats from an array of `VelocityVectorSphere`
"""
function get_coords_from_vel_array(vels::Array{VelocityVectorSphere})
    lats = [v.latd for v in vels]
    lons = [v.lond for v in vels]

    (lons, lats)
end

end # module
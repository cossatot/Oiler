module Utils

export predict_vels_from_poles, find_vel_cycles, diagonalize_matrices, 
    random_sample_vel_groups, build_Vc_from_vel_samples


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





function random_sample_vel(vel::VelocityVectorSphere)
    rnd = randn(2)
    ve = vel.ve + rnd[1] * vel.ee
    vn = vel.vn + rnd[2] * vel.en
    (ve, vn)
end


function random_sample_vels(vels::Array{VelocityVectorSphere}, n_samps::Int)
    rnd_ve_block = randn((length(vels), n_samps))
    rnd_vn_block = randn((length(vels), n_samps))

    ve_out = zeros(size(rnd_ve_block))
    vn_out = zeros(size(rnd_vn_block))

    for (row, vel) in enumerate(vels)
        ve_out[row,:] = vel.ve .+ rnd_ve_block[row,:] .* vel.ee
        vn_out[row,:] = vel.vn .+ rnd_ve_block[row,:] .* vel.en
    end
    (ve_out, vn_out)
end


function random_sample_vel_groups(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}},
n_samps::Int)
    vel_group_samps = Dict{Tuple{String,String},Dict{String,Array{Float64,2}}}()
    for (key, group) in vel_groups
        (ve, vn) = random_sample_vels(group, n_samps)
        vel_group_samps[key] = Dict("ve" => ve, "vn" => vn)
    end
    vel_group_samps
end


function Vc_trip_from_vals(ve::Float64, vn::Float64)
    [ve; vn; 0]
end


function build_Vc_from_vel_sample(vel_samp::Dict{String,Array{Float64,2}},
    ind::Int)

    reduce(vcat, [Vc_trip_from_vals(ve, vel_samp["vn"][i,ind])
    for (i, ve) in enumerate(vel_samp["ve"][:,ind])])

end


function build_Vc_from_vel_samples(vel_samps::Dict{Tuple{String,String},Dict{String,Array{Float64,2}}},
    vel_keys::Array{Tuple{String,String}}, ind::Int)

    reduce(vcat, [build_Vc_from_vel_sample(vel_samps[key], ind) for key in vel_keys])
end


flat(x, y = vcat(x...)) = x == y ? x : flat(y)

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


end # module
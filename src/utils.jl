#module Utils

#using ..BlockRotations

#include("./block_rotations.jl")
include("./io.jl")

#using Base.Iterators
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


function make_digraph_from_poles(poles::Union{Array{EulerPoleCart},Array{EulerPoleSphere}})
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


function build_constraint_matrix(cycle, vel_group_keys)
    constraint_mat = spzeros(3, length(vel_group_keys) * 3)

    for (v, inds) in cycle
        end_ind = 3 * inds["ind"]

        constraint_mat[3, end_ind] = inds["val"]
        constraint_mat[2, end_ind - 1] = inds["val"]
        constraint_mat[1, end_ind - 2] = inds["val"]
    end
    constraint_mat
end


function build_constraint_matrices(cycles, vel_group_keys)
    reduce(vcat, [build_constraint_matrix(cyc, vel_group_keys) for (i, cyc) in
    cycles])
end


function make_block_PvGb_from_vels(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}})

    vel_group_list = keys(vel_groups)

    big_PvGb = diagonalize_matrices([build_PvGb_from_vels(vel_groups[gr]) for gr
    in vel_group_list])
    
    return Dict("PvGb" => big_PvGb, 
                "keys" => collect(Tuple(keys(vel_groups))))
end


function make_block_inversion_matrices_from_vels(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}})

    out_dict = make_block_PvGb_from_vels(vel_groups)

    out_dict["Vc"] = reduce(vcat,
        [build_vel_column_from_vels(vel_groups[gr]) for gr in out_dict["keys"]])

    return out_dict
end

function make_block_inv_rhs(PvGb::SparseMatrixCSC{Float64,Int64},
Vc::Array{Float64,1}, constraint_rhs::Vector{Float64})
    rhs = [2 * PvGb' * Vc; constraint_rhs]
end


function make_block_inv_lhs_constraints(PvGb, cm)
    m, n = size(PvGb)
    p = size(cm)[1]

    lhs_term_1 = 2 * PvGb' * PvGb
    lhs = [lhs_term_1 cm'; cm spzeros(p, p)]
end


function set_up_block_inv_w_constraints(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}})

    vd = make_block_inversion_matrices_from_vels(vel_groups)
    cycles = find_vel_cycles(vd["keys"])

    cm = build_constraint_matrices(cycles, vd["keys"])
    p = size(cm)[1]

    lhs = make_block_inv_lhs_constraints(vd["PvGb"], cm)
    
    constraint_rhs = zeros(p)
    rhs = [2 * vd["PvGb"]' * vd["Vc"]; constraint_rhs]

    Dict("lhs" => lhs, "rhs" => rhs, "keys" => vd["keys"])
end


function solve_block_invs_from_vel_groups(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}})
    block_inv_setup = set_up_block_inv_w_constraints(vel_groups)

    kkt_soln = block_inv_setup["lhs"] \ block_inv_setup["rhs"]

    poles = Dict()
    for (i, (fix, mov)) in enumerate(block_inv_setup["keys"])
        poles[(fix, mov)] = EulerPoleCart(x = kkt_soln[i * 3 - 2], 
                                         y = kkt_soln[i * 3 - 1],
                                         z = kkt_soln[i * 3],
                                         fix = fix, mov = mov)
    end
    poles
end


function solve_for_block_poles_iterative(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}},
    n_iters::Int)
    # consider making .oiler_config file w/ defaults such as max_iters_per_chunk

    # left hand side
    vd = make_block_PvGb_from_vels(vel_groups)
    cycles = find_vel_cycles(vd["keys"])
    cm = build_constraint_matrices(cycles, vd["keys"])
    p = size(cm)[1]

    lhs = make_block_inv_lhs_constraints(vd["PvGb"], cm)
    lhs = lu(lhs)

    rand_vels = random_sample_vel_groups(vel_groups, n_iters)

    results = Dict("vels" => rand_vels, "poles" => Dict{Int,Array{EulerPoleCart,1}}())

    for iter in 1:n_iters
        Vc = build_Vc_from_vel_samples(rand_vels, vd["keys"], iter)
        rhs = make_block_inv_rhs(vd["PvGb"], Vc, zeros(p))

        kkt_soln = lhs \ rhs

        results["poles"][iter] = [EulerPoleCart(x = kkt_soln[i * 3 - 2],
                                                y = kkt_soln[i * 3 - 1],
                                                z = kkt_soln[i * 3],
                                                fix = fix, mov = mov)
               for (i, (fix, mov)) in enumerate(vd["keys"])]
    end
    results
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


flat(x,y = vcat(x...)) = x == y ? x : flat(y)

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


function get_pole_path(poles::Array{EulerPoleCart}, path::Array{String})

    steps = length(path) - 1
    pole_path = Array{EulerPoleCart,1}(undef, steps)


    for i in 1:steps
        place = path[i]
        next = path[i+1]

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

function get_path_euler_pole(poles::Array{EulerPoleCart,1}, fix::String, mov::String)

    vel_dg = make_digraph_from_poles(poles)
    vel_ug = make_ugraph_from_digraph(vel_dg)

    path = find_shortest_path(vel_ug, fix, mov)

    pole_path = get_pole_path(poles, path)

    final_pole = add_poles(pole_path)
end


#function flattenall(a::Array{Str})
#    while any(x->typeof(x)<:AbstractArray, a)
#        a = collect(Base.flatten(a))
#    end
#    return a
#end

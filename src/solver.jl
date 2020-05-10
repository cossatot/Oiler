module Solver

export make_block_PvGb_from_vels, solve_block_invs_from_vel_groups,
    solve_for_block_poles_iterative,
    make_block_inversion_matrices_from_vels
    
using ..Oiler
using ..Oiler: VelocityVectorSphere, PoleCart, PoleSphere, build_PvGb_from_vels,
    build_vel_column_from_vels, add_poles, pole_sphere_to_cart, find_vel_cycles,
    diagonalize_matrices, Fault

using Logging
using DataFrames
using SparseArrays
using LinearAlgebra

"""
    build_constraint_matrix(cycle, vel_group_keys)

Builds a matrix to add an equality constraint to the velocity solution. The
equality constraint is that a cycle of three graph-adjacent poles must sum
to zero, i.e. a_b + b_c + c_a = 0.
"""
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

"""
    build_constraint_matrices(cycle, vel_group_keys)

Builds a matrix to add equality constraints to the velocity solution. The
equality constraint is that a cycle of three graph-adjacent poles must sum
to zero, i.e. a_b + b_c + c_a = 0.  This function loops over all of the
velocity triangles (cycles) that are present in the pole adjacency graph.
"""
function build_constraint_matrices(cycles, vel_group_keys)
    reduce(vcat, [build_constraint_matrix(cyc, vel_group_keys) for (i, cyc) in
    cycles])
end


"""
    weight_from_error(error, zero_err_weight)

Returns either the inverse of the error, or `zero_err_weight` if the weight
is zero. The latter defaults to 1e20.
"""
function weight_from_error(error::Float64; zero_err_weight::Float64 = 1e-10)
    if error == 0.
        error = zero_err_weight
    end
    weight = error^-2
end
    

function build_weight_vector_from_vel(vel::VelocityVectorSphere)
    [weight_from_error(e) for e in (vel.ee, vel.en, vel.eu)]
end


function build_weight_vector_from_vels(vels::Array{VelocityVectorSphere})
    W = reduce(vcat, [build_weight_vector_from_vel(vel) for vel in vels])

    if all(W .== W[1])
        W = ones(size(W))
    end
    W
end


function
build_weight_vectors(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}})
    weight_groups = Dict{Tuple{String,String},Array{Float64,1}}()

    for (key, vels) in vel_groups
        weight_groups[key] = build_weight_vector_from_vels(vels)
    end
    weight_groups
end


function make_block_PvGb_from_vels(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}})
    v_keys = sort(collect(Tuple(keys(vel_groups))))

    big_PvGb = diagonalize_matrices([build_PvGb_from_vels(vel_groups[gr]) for gr
        in v_keys])

    weights = build_weight_vectors(vel_groups)
    weight_vec = reduce(vcat, [weights[key] for key in v_keys]) 
    
    return Dict("PvGb" => big_PvGb, 
                "keys" => v_keys,
                "weights" => weight_vec)
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


function add_fault_locking_to_PvGb(faults::Array{Fault}, 
    vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}},
    PvGb::SparseMatrixCSC{Float64,Int64})
    PvGb = Matrix(PvGb)

    @info "   calculating locking effects"
    locking_partials = Oiler.Elastic.calc_locking_effects(faults, vel_groups)
    @info "   done calculating locking effects"

    @info "   adding to PvGb"
    for (part_idx, partials) in locking_partials
        PvGb[part_idx[1], part_idx[2]] = PvGb[part_idx[1], part_idx[2]] + partials
    end
    @info "   done adding to PvGb"
    PvGb = sparse(PvGb)
    PvGb
end


function weight_inv_matrices(PvGb_in, Vc_in, weights)
    N = PvGb_in' * sparse(diagm(weights))
    Vc = N * Vc_in
    PvGb = N * PvGb_in
    (PvGb, Vc)
end


function
set_up_block_inv_w_constraints(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}};
    faults::Array = [], weighted::Bool = true)

    @info " making block inversion matrices"
    vd = make_block_inversion_matrices_from_vels(vel_groups)
    @info " done making block inversion matrices"
    @info " finding vel cycles"
    cycles = find_vel_cycles(vd["keys"])
    @info " done finding vel cycles"

    if length(faults) > 0
    @info " doing locking"
        PvGb = add_fault_locking_to_PvGb(faults, vel_groups, vd["PvGb"])
    @info " done doing locking"
    else
        PvGb = vd["PvGb"]
    end

    if weighted == true
        N = PvGb' * sparse(diagm(vd["weights"]))
        Vc = N * vd["Vc"]
        PvGb = N * PvGb
    else
        Vc = vd["Vc"]
    end

    if cycles == Dict{Any,Any}()
        lhs = PvGb
        rhs = Vc

    else
        cm = build_constraint_matrices(cycles, vd["keys"])
        p = size(cm)[1]

        # lhs = make_block_inv_lhs_constraints(PvGb, cm)
        lhs = [PvGb; 1e12 .* cm];
    
        constraint_rhs = zeros(p)
        # rhs = make_block_inv_rhs(PvGb, Vc, constraint_rhs)
        rhs = [Vc; constraint_rhs]
    end

    # rhs = [rhs; zeros(size(rhs,1))]
    # lhs = [lhs ; tikh_lambda * I]

    Dict("lhs" => lhs, "rhs" => rhs, "keys" => vd["keys"])
end


function solve_block_invs_from_vel_groups(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}};
    faults::Array = [], weighted::Bool = true, max_dense_size = 1000)

    @info "setting up matrices"
    block_inv_setup = set_up_block_inv_w_constraints(vel_groups; 
        faults = faults, weighted = weighted)
    @info "done with setup"
    
    lhs = block_inv_setup["lhs"]

    # if size(lhs, 1) < max_dense_size
    #    lhs = Matrix(lhs)
    # end

    @info "solving"
    kkt_soln = lhs \ block_inv_setup["rhs"]
    @info "done solving"

    poles = Dict()
    for (i, (fix, mov)) in enumerate(block_inv_setup["keys"])
        poles[(fix, mov)] = PoleCart(x = kkt_soln[i * 3 - 2], 
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

    results = Dict("vels" => rand_vels, "poles" => Dict{Int,Array{PoleCart,1}}())

    for iter in 1:n_iters
        Vc = build_Vc_from_vel_samples(rand_vels, vd["keys"], iter)
        rhs = make_block_inv_rhs(vd["PvGb"], Vc, zeros(p))

        kkt_soln = lhs \ rhs

        results["poles"][iter] = [PoleCart(x = kkt_soln[i * 3 - 2],
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


function Vc_triple_from_vals(ve::Float64, vn::Float64)
    [ve; vn; 0]
end


function build_Vc_from_vel_sample(vel_samp::Dict{String,Array{Float64,2}},
    ind::Int)

    reduce(vcat, [Vc_triple_from_vals(ve, vel_samp["vn"][i,ind])
    for (i, ve) in enumerate(vel_samp["ve"][:,ind])])

end


function build_Vc_from_vel_samples(vel_samps::Dict{Tuple{String,String},Dict{String,Array{Float64,2}}},
    vel_keys::Array{Tuple{String,String}}, ind::Int)

    reduce(vcat, [build_Vc_from_vel_sample(vel_samps[key], ind) for key in vel_keys])
end



end # module

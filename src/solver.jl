module Solver

export make_block_PvGb_from_vels, solve_block_invs_from_vel_groups,
    solve_for_block_poles_iterative,
    make_block_inversion_matrices_from_vels
    
    
using ..Oiler: VelocityVectorSphere, PoleCart, PoleSphere, build_PvGb_from_vels,
    build_vel_column_from_vels, add_poles, pole_sphere_to_cart, find_vel_cycles,
    diagonalize_matrices, random_sample_vel_groups, build_Vc_from_vel_samples


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
function weight_from_error(error::Float64; zero_err_weight::Float64 = 1e20)
    if error == 0.
        weight = zero_err_weight
    else
        weight = 1. / error
    end
    weight
end
    

function build_weight_vector_from_vel(vel::VelocityVectorSphere)
    [weight_from_error(e) for e in (vel.en, vel.ee, vel.ed)]
end


function build_weight_vector_from_vels(vels::Array{VelocityVectorSphere})
    reduce(vcat, [build_weight_vector_from_vel(vel) for vel in vels])
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

    vel_group_list = keys(vel_groups)

    big_PvGb = diagonalize_matrices([build_PvGb_from_vels(vel_groups[gr]) for gr
        in vel_group_list])

    v_keys = collect(Tuple(keys(vel_groups)))

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


function
set_up_block_inv_w_constraints(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}};
weighted::Bool = true, tikh_lambda::Float64 = 0.)

    vd = make_block_inversion_matrices_from_vels(vel_groups)
    cycles = find_vel_cycles(vd["keys"])
    
    if weighted == true
        PvGb = vd["weights"] .* vd["PvGb"]
        Vc = vd["weights"] .* vd["Vc"]
    else
        PvGb = vd["PvGb"]
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
    weighted::Bool = true, max_dense_size = 1000)

    block_inv_setup = set_up_block_inv_w_constraints(vel_groups; weighted =
        weighted)
    
    lhs = block_inv_setup["lhs"]

    # if size(lhs, 1) < max_dense_size
    #    lhs = Matrix(lhs)
    # end

    kkt_soln = lhs \ block_inv_setup["rhs"]

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

end # module
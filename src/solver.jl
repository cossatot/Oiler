module Solver

export make_block_PvGb_from_vels, solve_block_invs_from_vel_groups,
    solve_for_block_poles_iterative,
    make_block_inversion_matrices_from_vels

using ..Oiler
using ..Oiler: VelocityVectorSphere, PoleCart, PoleSphere, build_PvGb_from_vels,
    build_vel_column_from_vels, add_poles, pole_sphere_to_cart, find_vel_cycles,
    diagonalize_matrices, Fault

using Logging
using Statistics
using DataFrames
using SuiteSparse
using SparseArrays
using LinearAlgebra
using Distributions
import Base.Threads.@threads

LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
#SuiteSparse.UMFPACK.umf_ctrl[8] = 0

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
        constraint_mat[2, end_ind-1] = inds["val"]
        constraint_mat[1, end_ind-2] = inds["val"]
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
is zero. The latter defaults to 1e-10.

Note: this is applicable for uncorrelated (independent) errors.
"""
function weight_from_error(error::Float64; zero_err_weight::Float64=1e2)
    if error == 0.0
        error = zero_err_weight
    end
    weight = error^-2
end


function build_weight_vector_from_vel(vel::VelocityVectorSphere)
    [weight_from_error(e) for e in (vel.ee, vel.en, vel.eu)]
end


function build_weight_vector_from_vels(vels::Array{VelocityVectorSphere})
    W = reduce(vcat, [build_weight_vector_from_vel(vel) for vel in vels])
    W
end


function build_weight_vectors(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}})
    weight_groups = Dict{Tuple{String,String},Array{Float64,1}}()
    v_keys = sort(collect(Tuple(keys(vel_groups))))

    for (key, vels) in vel_groups
        weight_groups[key] = build_weight_vector_from_vels(vels)
    end
    weight_vec = reduce(vcat, [weight_groups[key] for key in v_keys])

    if all(weight_vec .== weight_vec[1])
        weight_vec = ones(size(weight_vec))
    end

    weight_vec
end


function var_cov_from_vel(vel::VelocityVectorSphere; zero_err_weight::Float64=1e2)

    if vel.ee == 0.0
        e_var = zero_err_weight
    else
        e_var = vel.ee^2
    end

    if vel.en == 0.0
        n_var = zero_err_weight
    else
        n_var = vel.en^2
    end

    if vel.eu == 0.0
        u_var = zero_err_weight
    else
        u_var = vel.eu^2
    end

    [e_var vel.cen 0.0; vel.cen n_var 0.0; 0.0 0.0 u_var]
end


function var_cov_from_tri(tri, zero_err_weight=1e-2)
    if tri.dip_slip_err == 0.0
        ds_var = zero_err_weight
    else
        ds_var = tri.dip_slip_err^2
    end

    if tri.strike_slip_err == 0.0
        ss_var = zero_err_weight
    else
        ss_var = tri.strike_slip_err^2
    end

    [ds_var tri.cds; tri.cds ss_var]
end


function build_var_cov_matrix_from_vels(vels::Array{VelocityVectorSphere})
    diagonalize_matrices([var_cov_from_vel(vel) for vel in vels])
end


function build_var_cov_weight_matrix(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}})
    v_keys = sort(collect(Tuple(keys(vel_groups))))
    weight_groups = Dict()

    for (key, vels) in vel_groups
        weight_groups[key] = build_var_cov_matrix_from_vels(vels)
    end
    diagonalize_matrices([weight_groups[key] for key in v_keys])
end


function make_block_PvGb_from_vels(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}})
    v_keys = sort(collect(Tuple(keys(vel_groups))))

    big_PvGb = diagonalize_matrices([build_PvGb_from_vels(vel_groups[gr]) for gr
                                     in
                                     v_keys])

    weight_vec = build_weight_vectors(vel_groups)

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
    PvGb::SparseMatrixCSC{Float64,Int64}; elastic_floor=1e-4, check_nans=false)

    @info "   calculating locking effects"
    @time locking_partials = Oiler.Elastic.calc_locking_effects(faults, vel_groups;
        elastic_floor=elastic_floor, check_nans=check_nans)
    #println(locking_partials)

    @info "   adding to PvGb"
    PvGb = Matrix(PvGb)
    @time for (part_idx, partials) in locking_partials
        PvGb[part_idx[1], part_idx[2]] = PvGb[part_idx[1], part_idx[2]] + partials
    end
    PvGb = sparse(PvGb)
    PvGb
end


function make_tri_regularization_matrix(tris, distance_weight)
    tri_eqns = []
    tri_adj_dict = Oiler.Tris.get_tri_adjacence_dict(tris)
    tri_enum = Dict([tri.name => i for (i, tri) in enumerate(tris)])

    for tri1 in tris
        tri1_name = tri1.name
        for tri2_name in tri_adj_dict[tri1.name]
            tri2 = tris[tri_enum[tri2_name]]
            t1_ind = 2 * tri_enum[tri1_name] - 1
            t2_ind = 2 * tri_enum[tri2_name] - 1

            const_mat = zeros(2, 2 * size(tris, 1))
            centroid_dist = Oiler.Tris.tri_centroid_distance(tri1, tri2)
            const_mat[1, t1_ind] = distance_weight / centroid_dist
            const_mat[2, t1_ind+1] = distance_weight / centroid_dist
            const_mat[1, t2_ind] = -distance_weight / centroid_dist
            const_mat[2, t2_ind+1] = -distance_weight / centroid_dist
            push!(tri_eqns, const_mat)
        end
    end
    vcat(tri_eqns...)
end


function make_tri_prior_matrices_old(tris)
    lhs = diagm(ones(2 * length(tris)))
    vels = ones(2 * length(tris))
    #weights = ones(2 * length(tris))

    for (i, tri) in enumerate(tris)
        ds_ind = 2 * i - 1
        ss_ind = 2 * i

        vels[ds_ind] = tri.dip_slip_rate
        #weights[ds_ind] = weight_from_error(tri.dip_slip_err^-1)
        vels[ss_ind] = tri.strike_slip_rate
        #weights[ss_ind] = weight_from_error(tri.strike_slip_err^-1)
    end
    weights = diagonalize_matrices(
        [var_cov_from_tri(tri) for tri in tris]
    )
    lhs, vels, weights
end


function make_tri_prior_matrices(tris)
    function is_default(tri)
        rates = [tri.dip_slip_rate, tri.dip_slip_err, tri.strike_slip_rate, tri.strike_slip_err]
        all(rates .== 0.0)
    end

    function get_tri_ind(tri, tris)
        for (i, t) in enumerate(tris)
            if t == tri
                return i
            end
        end
    end

    non_default_tris = filter(!is_default, tris)

    nt = length(tris)
    np = length(non_default_tris)
    lhs = zeros(np * 2, nt* 2)
    rhs = zeros(size(lhs)[1])
    var_covs = []

    for (i, tri) in enumerate(tris) #cols
        if !is_default(tri)
            ds_col_ind = 2 * i - 1
            ss_col_ind = 2 * i

            nd_tri_idx = get_tri_ind(tri, non_default_tris)
            ds_row_ind = 2 * nd_tri_idx - 1
            ss_row_ind = ds_row_ind + 1

            lhs[ds_row_ind, ds_col_ind] = 1.0
            lhs[ss_row_ind, ss_col_ind] = 1.0

            rhs[ds_row_ind] = tri.dip_slip_rate
            rhs[ss_row_ind] = tri.strike_slip_rate

            push!(var_covs, var_cov_from_tri(tri))
        end
    end
    weights = diagonalize_matrices(var_covs)
    lhs, rhs, weights
end



function make_tri_priors_from_rake_constraints(tri_rake_constraints, tris)
    n_tris = length(tris)
    n_constraints = length(tri_rake_constraints)
    
    lhs = zeros(2 * n_constraints, 2 * n_tris)
    vels = ones(2 * n_constraints)
    weights = Oiler.Utils.diagonalize_matrices(
        [var_cov_from_tri(trc) for trc in tri_rake_constraints]
        )
    
    for (i, trc) in enumerate(tri_rake_constraints)
        ss_row_ind = 2 * i
        ds_row_ind = ss_row_ind - 1

        tri_idx, tri = Oiler.IO.get_tri_from_tris(tris, trc.name)
        ss_col_ind = 2 * tri_idx
        ds_col_ind = ss_col_ind - 1

        lhs[ss_row_ind, ss_col_ind] += 1.0
        lhs[ds_row_ind, ds_col_ind] += 1.0

        vels[ss_row_ind] = trc.strike_slip_rate
        vels[ds_row_ind] = trc.dip_slip_rate
    end
    lhs, vels, weights
end


function add_tris_to_PvGb(tris, vel_groups, vd; priors=false, regularize=true,
        tri_rake_constraints=[], distance_weight::Float64=10.0, elastic_floor=1e-4)

    nt = length(tris)

    # Elastic locking effects on GNSS
    PvGb = vd["PvGb"]
    PvGb_n_cols = size(PvGb, 2)
    gnss_vels = get_gnss_vels(vel_groups)
    gnss_lons = [vel["vel"].lon for vel in gnss_vels]
    gnss_lats = [vel["vel"].lat for vel in gnss_vels]
    gnss_idxs = [vel["idx"] for vel in gnss_vels]

    @info "   calculating tri locking partials"
    @time tri_effects = Oiler.Elastic.calc_tri_effects(tris, gnss_lons, gnss_lats;
        elastic_floor=elastic_floor)
    @info "   done"
    #@info "    making tri effect matrices"
    tri_gnss_effects_matrix = zeros((size(PvGb)[1], size(tri_effects)[2]))

    #@info "    adding matrices together"
    for (i, vel_idx) in enumerate(gnss_idxs)
        i3 = i * 3
        gnss_row_idxs = i3-2:i3
        pvgb_row_idxs = vel_idx[1]

        tri_gnss_effects_matrix[pvgb_row_idxs, :] = tri_effects[gnss_row_idxs, :]
    end
    PvGb = hcat(PvGb, tri_gnss_effects_matrix)

    #@info "    doing regularization"
    # Priors and regularization
    tri_reg_matrix = zeros(0, nt * 2)
    #tri_reg_weights = zeros(0)
    tri_reg_weights = ones(0)
    tri_reg_vels = zeros(0)

    # @info "    making regularization matrices"
    if regularize == true
        distance_reg_matrix = make_tri_regularization_matrix(tris, distance_weight)
        tri_reg_matrix = vcat(tri_reg_matrix, distance_reg_matrix)
        tri_reg_weights = diagm(vcat(tri_reg_weights, ones(size(tri_reg_matrix, 1))))
        tri_reg_vels = vcat(tri_reg_vels, zeros(size(tri_reg_matrix, 1)))
    end

    #@info "    doing priors"
    if priors == true
        priors, prior_vels, prior_weights = make_tri_prior_matrices(tris)
        tri_reg_matrix = vcat(tri_reg_matrix, priors)
        tri_reg_vels = vcat(tri_reg_vels, prior_vels)
        tri_reg_weights = diagonalize_matrices([tri_reg_weights, prior_weights])
    else
        prior_weights = ones(size(tri_reg_weights))
    end

    if length(tri_rake_constraints) > 0
        rake_matrix, rake_vels, rake_weights = make_tri_priors_from_rake_constraints(
                                                    tri_rake_constraints, tris)
        #fake_rake_weights = ones(length(tri_rake_constraints))
        tri_reg_matrix = vcat(tri_reg_matrix, rake_matrix)
        #tri_reg_weights = vcat(tri_reg_weights, fake_rake_weights)
        tri_reg_weights = diagonalize_matrices([tri_reg_weights, rake_weights])
        tri_reg_vels = vcat(tri_reg_vels, rake_vels)
    end

    #@info "    finalizing"
    tri_reg_block = hcat(zeros(size(tri_reg_matrix, 1), PvGb_n_cols), tri_reg_matrix)
    PvGb = vcat(PvGb, tri_reg_block)

    #@info "    done"
    #vd["weights"] = vcat(vd["weights"], prior_weights)
    #vd["weights"] = vcat(vd["weights"], zeros(size(tri_reg_weights))) # maybe 1s?
    vd["tri_weights"] = tri_reg_weights
    vd["Vc"] = vcat(vd["Vc"], tri_reg_vels)
    vd["PvGb"] = PvGb

    vd
end


function weight_inv_matrices(PvGb_in, Vc_in, weights)
    if typeof(weights) == Vector{Float64}
        weights = sparse(diagm(weights))
    end
    N = PvGb_in' * weights
    Vc = N * Vc_in
    PvGb = N * PvGb_in
    (PvGb, Vc)
end


function add_equality_constraints_kkt(PvGb, Vc, cm)
    p = size(cm, 1)
    lhs = [2*PvGb'*PvGb cm'
        cm zeros(p, p)]
    rhs = [2 * PvGb' * Vc; zeros(p)]
    lhs, rhs
end


function add_equality_constraints_bi_objective(PvGb, Vc, cm)
    p = size(cm, 1)
    if p == 0
        lhs = PvGb
        rhs = Vc
    else
        lhs = [PvGb; cm .* 1e10]
        rhs = [Vc; zeros(p)]
    end

    (lhs, rhs)
end


function make_weighted_constrained_kkt_lls_matrices(PvGb, Vc, cm, weights; sparse_lhs::Bool=false)
    # W = sparse(diagm(1 ./ weights))
    W = weights
    p, q = size(cm)
    n = length(Vc)
    if sparse_lhs
        _zeros = spzeros
    else
        _zeros = zeros
    end
    lhs = [PvGb _zeros(n, p) W
        cm _zeros(p, p) _zeros(p, n)
        _zeros(q, q) cm' PvGb']

    rhs = [Vc; zeros(p); zeros(q)]

    if sparse_lhs
        return sparse(lhs), rhs
    else
        return lhs, rhs
    end
end


function make_weighted_constrained_kkt_lls_matrices_symmetric(PvGb, Vc, cm, weights; sparse_lhs::Bool=false)
    W = weights

    p, q = size(cm)
    n = length(Vc)
    if sparse_lhs
        _zeros = spzeros
    else
        _zeros = zeros
    end

    lhs = [_zeros(p, p) _zeros(p, n) cm
        _zeros(n, p) W PvGb
        cm' PvGb' _zeros(q, q)]

    rhs = [zeros(p); Vc; zeros(q)]

    if sparse_lhs
        return sparse(lhs), rhs
    else
        return lhs, rhs
    end
end


function make_weighted_constrained_lls_matrices(PvGb, Vc, cm, weights;
    sparse_lhs::Bool=false, method="kkt_sym")

    if method == "kkt"
        return make_weighted_constrained_kkt_lls_matrices(PvGb, Vc, cm, weights;
            sparse_lhs=sparse_lhs)
    elseif method == "kkt_sym"
        return make_weighted_constrained_kkt_lls_matrices_symmetric(PvGb, Vc, cm, weights;
            sparse_lhs=sparse_lhs)
    else
        throw(ArgumentError("contraint method not implemented"))
    end
end


function set_up_block_inv_no_constraints(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}};
        faults::Array=[], tris::Array=[], regularize_tris=true, tri_priors=false, 
        tri_rake_constraints=[],
    tri_distance_weight::Float64=10.0,
    elastic_floor=1e-4, check_nans=false, fill_nans=true, nan_fill_val=0.0)

    @info " making block inversion matrices"
    @time vd = make_block_inversion_matrices_from_vels(vel_groups)
    PvGb_size = size(vd["PvGb"])
    @info "  raw PvGb size: $PvGb_size"

    if check_nans == true
        @info " checking for NaNs"
        NaNs = false
        if any(isnan, vd["PvGb"])
            @warn "NaNs in PvGb"
            NaNs = true
            if fill_nans == true
                for i = eachindex(vd["PvGb"])
                    if isnan(vd["PvGb"][i])
                        vd["PvGb"][i] = nan_fill_val
                    end
                end
            end
        end
        if any(isnan, vd["Vc"])
            @warn "NaNs in Vc"
            NaNs = true
            if fill_nans == true
                for i = eachindex(vd["Vc"])
                    if isnan(vd["Vc"][i])
                        vd["Vc"][i] = nan_fill_val
                    end
                end
            end
        end
        if ~NaNs
            @info "No NaNs in matrices"
        end
    end


    if length(faults) > 0
        @info " doing locking"
        vd["PvGb"] = add_fault_locking_to_PvGb(faults, vel_groups, vd["PvGb"];
            elastic_floor=elastic_floor, check_nans=check_nans)
        @info " done doing locking"

    end

    if check_nans == true
        @info " checking for NaNs"
        NaNs = false
        if any(isnan, vd["PvGb"])
            @warn "NaNs in PvGb"
            NaNs = true
        end
        if any(isnan, vd["Vc"])
            @warn "NaNs in Vc"
            NaNs = true
        end
        if ~NaNs
            @info "No NaNs in matrices"
        end
    end

    if length(tris) > 0
        @info " doing tris"
        @time vd = add_tris_to_PvGb(tris, vel_groups, vd;
            priors=tri_priors, regularize=regularize_tris,
            tri_rake_constraints=tri_rake_constraints,
            distance_weight=tri_distance_weight,
            elastic_floor=elastic_floor)
        @info " done doing tris"
    end

    if check_nans == true
        @info " checking for NaNs"
        NaNs = false
        if any(isnan, vd["PvGb"])
            @warn "NaNs in PvGb"
            NaNs = true
            if fill_nans == true
                for i = eachindex(vd["PvGb"])
                    if isnan(vd["PvGb"][i])
                        vd["PvGb"][i] = nan_fill_val
                    end
                end
            end
            if any(isnan, vd["PvGb"])
                @warn "Couldn't replace NaNs in PvGb"
            end
        end
        if any(isnan, vd["Vc"])
            @warn "NaNs in Vc"
            NaNs = true
            if fill_nans == true
                for i = eachindex(vd["Vc"])
                    if isnan(vd["Vc"][i])
                        vd["Vc"][i] = nan_fill_val
                    end
                end
            end
            if any(isnan, vd["Vc"])
                @warn "Couldn't replace NaNs in Vc"
            end
        end
        if ~NaNs
            @info "No NaNs in matrices"
        end
    end

    @info " done making block inversion matrices"
    vd
end


function set_up_block_inv_w_constraints(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}};
        basic_matrices::Dict=Dict(), faults::Array=[], tris::Array=[], tri_rake_constraints=[],
    weighted::Bool=true,
    elastic_floor=1e-4,
    regularize_tris::Bool=true, tri_priors::Bool=false,
    tri_distance_weight::Float64=10.0, sparse_lhs::Bool=false,
    constraint_method="kkt_sym", check_nans=false)

    if basic_matrices == Dict()
        basic_matrices = set_up_block_inv_no_constraints(vel_groups;
            faults=faults, tris=tris, regularize_tris=regularize_tris,
            tri_priors=tri_priors, tri_rake_constraints=tri_rake_constraints,
            tri_distance_weight=tri_distance_weight,
            elastic_floor=elastic_floor, check_nans=check_nans)
    end

    PvGb = basic_matrices["PvGb"]
    Vc = basic_matrices["Vc"]
    weights = basic_matrices["weights"]

    @info " finding vel cycles"
    @time cycles = find_vel_cycles(basic_matrices["keys"])
    @info " done finding vel cycles"

    if cycles == Dict{Any,Any}()
        # lhs = PvGb
        # rhs = Vc
        cm = []
    else
        cm = build_constraint_matrices(cycles, basic_matrices["keys"])
        if length(tris) > 0
            cm = hcat(cm, zeros(size(cm)[1], length(tris) * 2))
            cm = vcat(cm, zeros(size(PvGb, 1) - size(cm, 1), size(PvGb, 2)))
        end
        cm = Oiler.Utils.sort_sparse_matrix(cm)
    end

    if weighted == true
        if cycles == Dict{Any,Any}()
            @info " Using weighted, non-constrained LLS"
            #@warn "does not currently use velocity covariances"
            weights = build_var_cov_weight_matrix(vel_groups)
            lhs, rhs = weight_inv_matrices(PvGb, Vc, weights)
        else
            @info " Using weighted, constrained LLS"
            cm = cm[Oiler.Utils.lin_indep_rows(cm), :]
            weights = build_var_cov_weight_matrix(vel_groups)
            basic_matrices["var_cov_matrix"] = weights

            if length(tris) > 0
                weights = diagonalize_matrices([weights,
                    #(basic_matrices["tri_weights"] .^ 2)])
                    (basic_matrices["tri_weights"])])
            end

            lhs, rhs = make_weighted_constrained_lls_matrices(PvGb, Vc, cm,
                weights; sparse_lhs=sparse_lhs, method=constraint_method)
        end
    else
        if length(cm) == 0
            @info " Using non-weighted non-constrained LLS"
            lhs, rhs = PvGb, Vc
        else
            cm = cm[Oiler.Utils.lin_indep_rows(cm), :]
            @info " Using non-weighted, constrained LLS"
            lhs, rhs = add_equality_constraints_kkt(PvGb, Vc, cm)
        end
    end

    basic_matrices["lhs"] = lhs
    basic_matrices["rhs"] = rhs
    basic_matrices["cm"] = cm

    basic_matrices
end

"""
    solve_block_invs_from_vel_groups

Main solution function in Oiler. This function takes the `vel_groups` dictionary
that holds all of the velocity information as well as any faults or tris, and
solves for Euler poles, fault/tri slip rates, and uncertainty.

# Arguments
- `vel_groups`: Dictionary with lists of `VelocityVectorSphere` organized by the
    pole that the velocity is referenced to.
- `faults`: (Optional) list of `Fault`s that are used in the inversion. Each
    `Fault` is used to calculate geodetic locking effects and will provide an
    additional pole between the hanging wall and footwall block, which can
    change a model from an unconstrained to a constrained least squares model. 
- `tris`: (Optional) list of `Tri`s that are used to simulate large, irregular
    fault surfaces that may slip at a rate different than the block convergence
    rate. Tris are used to model geodetic locking but (as they generally do not
    slip at the block convergence rate) do not provide an additional pole.
- `weighted`: Boolean flag to indicate whether to incorporate observation
    weights (i.e. standard deviations of the observed velocity components) into
    the solution.
- `regularize`: (*semi-deprecated*) Boolean flag to indicate whether to incorporate L2
    regularization into the solution. It is (probably) only applicable for
    models without equality constraints on the poles, i.e. without faults,
    because it does not seem to be compatible with the Constrained Least Squares
    solution.
- `l2_lambda`: (*semi-deprecated*) Value to use for L2 regularization. Defaults
    to 100.
- `check_closures`: Boolean flag to indicate whether to check that plate
    circuits close. Defaults to `true`.
- `sparse_lhs`: Boolean flag to indicate whether the design matrix (left-hand
    side or LHS in the solution) should be sparse or dense. A sparse matrix will
    probably take less RAM and use the sparse solver from SuiteSparse, which may
    perform differently than the dense solver. Defaults to `false`.
- `predict_vels`: Boolean flag to indicate whether to predict GNSS velocities
    and fault slip rates and deliver them in the results. Defaults to `false`.
- `constraint_method`: Method used to deal with pole circuits (equality
    constraints) in the solution. Currently the only method that is known to
    work is making a system that uses Karush-Kuhn-Tucker conditions.  Defaults
    to `"kkt"`.
- `pred_se`: Boolean flag whether to calculate standard errors for the poles and
    slip rates. This can be more computationally expensive than the regular
    solution, and may not need to be run regularly. Defaults to `false`.

"""
function solve_block_invs_from_vel_groups(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}};
    faults::Array=[], tris=[], weighted::Bool=true,
    regularize_tris::Bool=true,
    tri_priors::Bool=false, 
    tri_rake_constraints=[],
    elastic_floor=1e-4,
    tri_distance_weight::Float64=10.0, check_closures::Bool=true,
    sparse_lhs::Bool=false, predict_vels::Bool=false, constraint_method="kkt_sym",
    pred_se::Bool=false, se_iters=1000, factorization="lu", check_nans=false)

    @info "setting up unconstrained matrices"
    @time block_inv_setup = set_up_block_inv_no_constraints(vel_groups;
        faults=faults,
        regularize_tris=regularize_tris,
        tri_priors=tri_priors,
        tri_distance_weight=tri_distance_weight,
        tri_rake_constraints=tri_rake_constraints,
        elastic_floor=elastic_floor,
        tris=tris, check_nans=check_nans)
    @info "making constrained matrices"
    @time block_inv_setup = set_up_block_inv_w_constraints(vel_groups;
        basic_matrices=block_inv_setup,
        weighted=weighted,
        regularize_tris=regularize_tris,
        tri_priors=tri_priors,
        tris=tris,
        tri_distance_weight=tri_distance_weight,
        tri_rake_constraints=tri_rake_constraints,
        elastic_floor=elastic_floor,
        sparse_lhs=sparse_lhs,
        constraint_method=constraint_method,
        check_nans=check_nans)

    lhs = block_inv_setup["lhs"]
    lhs_size = size(lhs)
    rhs = block_inv_setup["rhs"]
    rhs_size = size(block_inv_setup["rhs"])
    @info "LHS block size $lhs_size"
    @info "RHS size $rhs_size"

    @info uppercase(factorization) * " factorization"
    if factorization == "lu"
        fact = lu
    elseif factorization == "svd"
        fact = svd
        # ensure dense LHS
        lhs = Matrix(lhs)
    elseif factorization == "qr"
        fact = qr
    end
    @time lhs_fact = fact(lhs)
    @info "  Done with factorization"

    @info "solving"
    @time full_soln = lhs_fact \ rhs

    nk = length(full_soln)
    np = length(block_inv_setup["keys"])
    nt = length(tris)

    soln_length = np * 3 + nt * 2
    if constraint_method == "kkt_sym"
        soln_idx = nk-soln_length+1:nk
    else
        soln_idx = 1:soln_length
    end

    tri_soln_idx = (np*3+1):(np*3)+nt*2
    soln = full_soln[soln_idx]
    tri_soln = soln[tri_soln_idx]

    results = Dict()
    results["poles"] = get_poles_from_soln(block_inv_setup["keys"], soln)
    results["tri_slip_rates"] = get_tri_rates(tri_soln, tris)

    results["stats_info"] = Dict{Any,Any}(
        "RMSE" => Oiler.ResultsAnalysis.calc_RMSE_from_G(block_inv_setup, results)
    )
    RMSE_string = "RMSE: " * string(results["stats_info"]["RMSE"])
    @info RMSE_string
    #results["stats_info"]["n_obs"], results["stats_info"]["n_params"] = size(
    #    block_inv_setup["PvGb"]
    #)
    
    results["stats_info"]["n_obs"] = length(Oiler.Utils.get_gnss_vels(vel_groups)) + 
                                     length(Oiler.Utils.get_geol_slip_rate_vels(vel_groups))
    results["stats_info"]["n_params"] = 3 * Oiler.Utils.get_n_blocks_from_vel_groups(vel_groups) + nt
    


    if pred_se == true
        if weighted == true
            se_weights = block_inv_setup["weights"]
        else
            se_weights = ones(size(block_inv_setup["weights"]))
        end
        @info "Estimating solution uncertainties"
        @time pole_var = get_soln_covariance_matrix(block_inv_setup,
            lhs_fact,
            results,
            soln_idx,
            constraint_method;
            n_iters=se_iters,
            weighted=weighted
        )
        get_pole_uncertainties!(results["poles"], pole_var,
            block_inv_setup["keys"])
        if nt > 0
            get_tri_uncertainties!(pole_var, tris, results["tri_slip_rates"],
                tri_soln_idx)
        end
    end

    if predict_vels == true
        results["predicted_vels"] = Oiler.ResultsAnalysis.predict_model_velocities(
            vel_groups,
            block_inv_setup, results["poles"];
            tri_results=results["tri_slip_rates"])

        results["predicted_slip_rates"] = Oiler.ResultsAnalysis.predict_slip_rates(
            faults, results["poles"])
    end

    if check_closures == true
        closures = Oiler.Utils.check_vel_closures(results["poles"]; tol=1e-4)
    end

    results
end


function get_poles_from_soln(vel_group_keys, soln)
    poles = Dict()
    for (i, (fix, mov)) in enumerate(vel_group_keys)
        poles[(fix, mov)] = PoleCart(x=soln[i*3-2],
            y=soln[i*3-1],
            z=soln[i*3],
            fix=fix,
            mov=mov)
    end
    poles
end

function diag_dot(A, B)
    out = dot.(eachrow(A), eachcol(B))
end


function make_OLS_cov(PvGb)
    @info "OLS covariances"
    Matrix(PvGb * PvGb') \ I
end


function make_WLS_cov(PvGb, weights)
    @info "WLS covariances"
    Matrix(PvGb' * sparse(diagm(weights)) * PvGb) \ I
end


function make_CLS_cov(PvGb, cm)
    @info "CLS covariances"
    @warn "unclear if CLS covariance calcs are correct"
    n, p = size(PvGb)
    nc = size(cm, 1)

    J = [Matrix(PvGb); Matrix(cm)]
    J_inv = pinv(J)

    cov = J_inv * [I(n) zeros(n, nc); zeros(nc, n) zeros(nc, nc)] * J_inv'
end


function make_CWLS_cov_iter(lhs, var_cov_matrix, Vc, cm, n_pole_vars, soln_idx,
    n_iters)
    @info "CWLS covariances through Monte Carlo techniques"
    stoch_poles = make_stoch_poles(lhs, var_cov_matrix, Vc, cm, n_pole_vars,
        soln_idx, n_iters, constraint_method)
    var_cov = cov(stoch_poles; dims=1)
end


function make_stoch_poles(lhs, var_cov_matrix, Vc, cm, n_pole_vars, soln_idx,
    n_iters, constraint_method)
    # this one is for CWLS

    p, q = size(cm)

    # rng = MersenneTwister(69) # to be replaced via config later
    # rand_vel_noise = randn((length(Vc), n_iters))

    # this works but is a bit RAM intensive
    vel_stds = rand(Distributions.MvNormal(Matrix(var_cov_matrix)), n_iters)

    n_remaining = length(Vc) - size(var_cov_matrix, 1)
    zeros_remaining = zeros(n_remaining)

    stoch_poles = zeros((n_iters, n_pole_vars)) # we want each soln to be a row

    @threads for i = 1:n_iters
        # Vc_stochastic = Vc + vel_stds .* rand_vel_noise[:,i]
        Vc_stochastic = Vc + vcat(vel_stds[:, i], zeros_remaining)
        if constraint_method == "kkt_sym"
            rhs = [zeros(p); Vc_stochastic; zeros(q)]
        else
            rhs = [Vc_stochastic; zeros(p); zeros(q)]
        end

        full_soln = lhs \ rhs
        soln = full_soln[soln_idx]
        stoch_poles[i, :] = soln'
    end
    stoch_poles
end


function get_soln_covariance_matrix(block_matrices, lhs_fact, results, soln_idx,
    constraint_method; n_iters=1000,
    weighted=true, save_stoch_poles=true)
    PvGb = block_matrices["PvGb"]
    cm = block_matrices["cm"]
    y_obs = block_matrices["Vc"]
    n, p = size(PvGb)

    if weighted == true
        weights = block_matrices["weights"]
    else
        weights = ones(length(block_matrices["weights"]))
    end

    if all(x -> x == 1.0, weights) & (length(cm) == 0)
        var_cov = make_OLS_cov(PvGb)
    elseif any(x -> x != 1.0, weights) & (length(cm) == 0)
        var_cov = make_WLS_cov(PvGb, weights)
    else
        stoch_poles = make_stoch_poles(lhs_fact,
            block_matrices["var_cov_matrix"], y_obs, cm,
            p, soln_idx, n_iters, constraint_method)
        var_cov = cov(stoch_poles; dims=1)
        if save_stoch_poles == true
            block_matrices["stoch_poles"] = stoch_poles
        end
    end

    var = results["stats_info"]["RMSE"]^2 * var_cov
    standard_error_vec = sqrt.(diag(var))

    SE_string = "mean standard error: " * string(
        sum(standard_error_vec) / length(standard_error_vec))
    @info SE_string

    var
end


function get_pole_uncertainties!(poles, pole_var, vel_group_keys)
    for (i, (fix, mov)) in enumerate(vel_group_keys)
        poles[(fix, mov)] = PoleCart(poles[(fix, mov)];
            ex=sqrt(pole_var[i*3-2, i*3-2]),
            ey=sqrt(pole_var[i*3-1, i*3-1]),
            ez=sqrt(pole_var[i*3, i*3]),
            cxy=pole_var[i*3-2, i*3-1],
            cxz=pole_var[i*3-2, i*3],
            cyz=pole_var[i*3-1, i*3])
    end
end


function get_tri_rates(tri_soln, tris)
    tri_results = Dict()

    if length(tris) > 0

        for (i, tri) in enumerate(tris)
            ds_ind = i * 2 - 1
            ss_ind = i * 2
            tri_results[tri.name] = Dict()
            tri_results[tri.name]["dip_slip"] = tri_soln[ds_ind]
            tri_results[tri.name]["strike_slip"] = tri_soln[ss_ind]
        end
    end
    tri_results
end


function get_tri_uncertainties!(pole_var, tris, tri_results, tri_soln_idx)
    tri_var = pole_var[tri_soln_idx, tri_soln_idx]
    if length(tris) > 0
        for (i, tri) in enumerate(tris)
            ds_ind = i * 2 - 1
            ss_ind = i * 2
            tri_results[tri.name]["dip_slip_err"] = sqrt(tri_var[ds_ind, ds_ind])
            tri_results[tri.name]["strike_slip_err"] = sqrt(tri_var[ss_ind, ss_ind])
            tri_results[tri.name]["dsc"] = tri_var[ds_ind, ss_ind]
        end
    end
end



function calc_forward_velocities(forward_vels, fix; poles, stoch_poles=[], faults=[],
    tris=[], tri_soln=[], elastic_floor=0.0)
    # this should include locking (or not) for all vels

    forward_vel_groups = group_vels_by_fix_mov(forward_vels)

    fault_vels = Oiler.IO.make_vels_from_faults(faults)

    if length(faults) > 0
        vel_groups = group_vels_by_fix_mov(vcat(forward_vels, fault_vels))
    else
        vel_groups = forward_vel_groups
    end

    block_mats = make_block_PvGb_from_vels(vel_groups)

    results_poles = Dict(k => Oiler.Utils.get_path_euler_pole(poles, k[1], k[2])
                         for k in block_mats["keys"])

    if length(faults) > 0
        block_mats["PvGb"] = add_fault_locking_to_PvGb(faults, vel_groups,
            block_mats["PvGb"]; elastic_floor=elastic_floor)
    end

    if length(tris) > 0
        @warn "can't calculate forward vels including tris yet!"
        #block_mats = add_tris_to_PvGb(tris, forward_vel_groups, block_mats;
        #    regularize=false, priors=false, elastic_floor=elastic_floor)
    end

    # pred results dict like real results dict
    pred_results = Dict("predicted_vels" => Oiler.ResultsAnalysis.predict_model_velocities(vel_groups,
        block_mats, results_poles))

    if length(stoch_poles) > 0

    end

    pred_vels_at_sites_ = Oiler.ResultsAnalysis.get_gnss_results(pred_results, vel_groups)

    pred_vels_at_sites = pred_vels_at_sites_[!, [:lon, :lat, :pred_ve, :pred_vn]]

    return pred_vels_at_sites
end


function solve_for_block_poles_iterative(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}},
    n_iters::Int)

    @warn "not updated in a while and perhaps disfunctional"
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

    for iter = 1:n_iters
        Vc = build_Vc_from_vel_samples(rand_vels, vd["keys"], iter)
        rhs = make_block_inv_rhs(vd["PvGb"], Vc, zeros(p))

        kkt_soln = lhs \ rhs

        results["poles"][iter] = [PoleCart(x=kkt_soln[i*3-2],
            y=kkt_soln[i*3-1],
            z=kkt_soln[i*3],
            fix=fix, mov=mov)
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


function random_sample_vels(vels::Array{VelocityVectorSphere}, n_samps::Int; only_data_vels=true)
    rnd_ve_block = randn((length(vels), n_samps))
    rnd_vn_block = randn((length(vels), n_samps))

    ve_out = zeros(size(rnd_ve_block))
    vn_out = zeros(size(rnd_vn_block))

    not_data_vel_types = ["Boundary", "fault", "tri"]

    for (row, vel) in enumerate(vels)
        if only_data_vels == true && vel.vel_type in not_data_vel_types
            ve_out[row, :] = vel.ve
            vn_out[row, :] = vel.vn
        else
            ve_out[row, :] = vel.ve .+ rnd_ve_block[row, :] .* vel.ee
            vn_out[row, :] = vel.vn .+ rnd_ve_block[row, :] .* vel.en
        end
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

    reduce(vcat, [Vc_triple_from_vals(ve, vel_samp["vn"][i, ind])
                  for (i, ve) in enumerate(vel_samp["ve"][:, ind])])

end


function build_Vc_from_vel_samples(vel_samps::Dict{Tuple{String,String},Dict{String,Array{Float64,2}}},
    vel_keys::Array{Tuple{String,String}}, ind::Int)

    reduce(vcat, [build_Vc_from_vel_sample(vel_samps[key], ind) for key in vel_keys])
end



end # module

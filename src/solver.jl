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
import Base.Threads.@threads

LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
SuiteSparse.UMFPACK.umf_ctrl[8] = 0

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
is zero. The latter defaults to 1e-10.
"""
function weight_from_error(error::Float64; zero_err_weight::Float64=1e2)
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

    if all(weight_vec .== weight_vec[1])
        weight_vec = ones(size(weight_vec))
    end
    
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

    @info "   calculating locking effects"
    @time locking_partials = Oiler.Elastic.calc_locking_effects(faults, vel_groups)
    
    @info "   adding to PvGb"
    PvGb = Matrix(PvGb)
    @time for (part_idx, partials) in locking_partials
        PvGb[part_idx[1], part_idx[2]] = PvGb[part_idx[1], part_idx[2]] + partials
    end
    PvGb = sparse(PvGb)
    PvGb
end


function make_tri_regularization_matrix(tris; distance_weight=1e3)
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
            const_mat[1,t1_ind] = distance_weight / centroid_dist
            const_mat[2,t1_ind + 1] = distance_weight / centroid_dist
            const_mat[1,t2_ind] = -distance_weight / centroid_dist
            const_mat[2,t2_ind] = -distance_weight / centroid_dist
            push!(tri_eqns, const_mat)
        end
    end
    vcat(tri_eqns...)
end


function add_tris_to_PvGb(tris, vel_groups, vd)
    PvGb = vd["PvGb"]
    gnss_vels = get_gnss_vels(vel_groups)
    gnss_lons = [vel["vel"].lon for vel in gnss_vels]
    gnss_lats = [vel["vel"].lat for vel in gnss_vels]
    gnss_idxs = [vel["idx"] for vel in gnss_vels]
    
    tri_reg_matrix = make_tri_regularization_matrix(tris)
    
    tri_effects = Oiler.Elastic.calc_tri_effects(tris, gnss_lons, gnss_lats)
    tri_eqn_matrix = zeros((size(PvGb)[1], size(tri_effects)[2]))

    for (i, vel_idx) in enumerate(gnss_idxs)
        i3 = i * 3
        gnss_row_idxs = i3 - 2:i3
        pvgb_row_idxs = vel_idx[1]

        tri_eqn_matrix[pvgb_row_idxs, :] = tri_effects[gnss_row_idxs, :]
    end

    tri_reg_block = hcat(zeros(size(tri_reg_matrix, 1), size(PvGb, 2)), tri_reg_matrix)
    PvGb = hcat(PvGb, tri_eqn_matrix)
    PvGb = vcat(PvGb, tri_reg_block)
    vd["weights"] = vcat(vd["weights"], ones(size(tri_reg_matrix, 1)))
    vd["Vc"] = vcat(vd["Vc"], zeros(size(tri_reg_matrix, 1)))
    vd["PvGb"] = PvGb
    vd
end


function weight_inv_matrices(PvGb_in, Vc_in, weights)
    N = PvGb_in' * sparse(diagm(weights))
    Vc = N * Vc_in
    PvGb = N * PvGb_in
    (PvGb, Vc)
end


function add_equality_constraints_kkt(PvGb, Vc, cm)
    p = size(cm, 1)
    lhs = [2 * PvGb' * PvGb  cm'; 
           cm                zeros(p, p)]
    rhs = [2 * PvGb' * Vc;   zeros(p)]
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
    W = sparse(diagm(1 ./ weights))
    # W = sparse(diagm(weights))
    p, q = size(cm)
    n = length(Vc)
    if sparse_lhs
        _zeros = spzeros
    else
        _zeros = zeros
    end
    lhs = [PvGb         _zeros(n, p) W;
           cm           _zeros(p, p) _zeros(p, n);
           _zeros(q, q) cm'          PvGb']

    rhs = [Vc; zeros(p); zeros(q)]

    if sparse_lhs
        return sparse(lhs), rhs
    else
        return lhs, rhs
    end
end


function make_weighted_constrained_kkt_lls_matrices_symmetric(PvGb, Vc, cm, weights; sparse_lhs::Bool=false)
    W = sparse(diagm(1 ./ weights))
    # W = sparse(diagm(weights))
    p, q = size(cm)
    n = length(Vc)
    if sparse_lhs
        _zeros = spzeros
    else
        _zeros = zeros
    end
    # lhs = [PvGb         _zeros(n, p) W;
    #       cm           _zeros(p, p) _zeros(p, n);
    #       _zeros(q, q) cm'          PvGb']

    lhs = [_zeros(p, p) _zeros(p, n) cm;
           _zeros(n, p)  W           PvGb;
           cm'           PvGb'       _zeros(q, q)]

    rhs = [zeros(p); Vc; zeros(q)]

    # rhs = [Vc; zeros(p); zeros(q)]

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
        # return make_weighted_constrained_kkt_lls_matrices_symmetric(PvGb, Vc, cm, weights;
            sparse_lhs=sparse_lhs)
    elseif method == "kkt_sym"
        return make_weighted_constrained_kkt_lls_matrices_symmetric(PvGb, Vc, cm, weights;
            sparse_lhs=sparse_lhs)
    else
        throw(ArgumentError("contraint method not implemented"))
    end
end


function set_up_block_inv_no_constraints(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}};
    faults::Array=[], tris::Array=[])

    @info " making block inversion matrices"
    @time vd = make_block_inversion_matrices_from_vels(vel_groups)
    PvGb_size = size(vd["PvGb"])
    @info "  raw PvGb size: $PvGb_size"


    if length(faults) > 0
        @info " doing locking"
        vd["PvGb"] = add_fault_locking_to_PvGb(faults, vel_groups, vd["PvGb"])
        @info " done doing locking"

    end

    if length(tris) > 0
        @info " doing tris"
        @time vd = add_tris_to_PvGb(tris, vel_groups, vd)
        @info " done doing tris"
    end
    
    @info " done making block inversion matrices"
    vd
end


function set_up_block_inv_w_constraints(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}};
    basic_matrices::Dict=Dict(), faults::Array=[], tris::Array=[], weighted::Bool=true, 
    regularize::Bool=false, l2_lambda::Float64=100.0, sparse_lhs::Bool=false,
    constraint_method="kkt")

    if basic_matrices == Dict()
        basic_matrices = set_up_block_inv_no_constraints(vel_groups, faults=faults, tris=tris)
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
    end

    if regularize == true
        @info "regularizing: l2_lambda = $l2_lambda"
        n_vars = size(lhs, 2)
        reg_matrix = sparse(1:n_vars, 1:n_vars, ones(n_vars)) .* l2_lambda

        rhs_reg = zeros(n_vars)

        lhs = [lhs; reg_matrix]
        rhs = [rhs; rhs_reg]
    end
    
    if weighted == true
        if cycles == Dict{Any,Any}()
            @info " Using weighted, non-constrained LLS"
            lhs, rhs = weight_inv_matrices(PvGb, Vc, weights)
        else
            @info " Using weighted, constrained LLS"
            cm = cm[Oiler.Utils.lin_indep_rows(cm), :]
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
    faults::Array=[], tris=[], weighted::Bool=true, regularize::Bool=false,
    l2_lambda::Float64=100.0, check_closures::Bool=true, sparse_lhs::Bool=false,
    predict_vels::Bool=false, constraint_method="kkt", pred_se::Bool=false,
    factorization="lu")

    @info "setting up unconstrained matrices"
    @time block_inv_setup = set_up_block_inv_no_constraints(vel_groups;
            faults=faults,
            tris=tris)
    @info "making constrained matrices"
    @time block_inv_setup = set_up_block_inv_w_constraints(vel_groups;
            basic_matrices=block_inv_setup,
            weighted=weighted,
            tris=tris,
            regularize=regularize, l2_lambda=l2_lambda,
            sparse_lhs=sparse_lhs,
            constraint_method=constraint_method)
    
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
    end
    @time lhs_fact = fact(lhs)
    @info "  Done with factorization"

    @info "solving"
    @time full_soln = lhs_fact \ rhs
    
    nk = length(full_soln)
    np = length(block_inv_setup["keys"])
    nt = length(tris)
    
    # println("nk: $nk")
    # println("np: $np")
    # println("nt: $nt")

    soln_length = np * 3 + nt * 2
    if constraint_method == "kkt_sym"
        soln_idx = nk - soln_length + 1:nk
    else
        soln_idx = 1:soln_length
    end

    tri_soln_idx = (np * 3 + 1):(np * 3) + nt * 2
    soln = full_soln[soln_idx]
    # println("soln_idx: ", soln_idx)
    # println("tri_soln_idx: ", tri_soln_idx)

    results = Dict()
    results["poles"] = get_poles_from_soln(block_inv_setup["keys"], soln)
    tri_soln = soln[tri_soln_idx]
    results["tri_slip_rates"] = get_tri_rates(tri_soln, tris)

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
                                                    weighted,
                                                    )
        get_pole_uncertainties!(results["poles"], pole_var, 
                                block_inv_setup["keys"])
        if nt > 0
            get_tri_uncertainties!(pole_var, tris, results["tri_slip_rates"], 
                                   tri_soln_idx)
        end
    end

    if predict_vels == true
        results["predicted_vels"] = predict_model_velocities(vel_groups, 
            block_inv_setup, results["poles"]; 
            tri_results=results["tri_slip_rates"])
        
        results["predicted_slip_rates"] = predict_slip_rates(faults, 
                                                             results["poles"])
    end

    if check_closures == true
        closures = Oiler.Utils.check_vel_closures(results["poles"])
    end

    results
end


function get_poles_from_soln(vel_group_keys, soln)
    poles = Dict()
    for (i, (fix, mov)) in enumerate(vel_group_keys)
        poles[(fix, mov)] = PoleCart(x=soln[i * 3 - 2], 
                                     y=soln[i * 3 - 1],
                                     z=soln[i * 3],
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


function make_CWLS_cov_iter(lhs, weights, Vc, cm, n_pole_vars, soln_idx, n_iters)
    @info "CWLS covariances through Monte Carlo techniques"
    p, q = size(cm)
    zvec = [zeros(p); zeros(q)]

    vel_stds = sqrt.(1 ./ weights)

    # rng = MersenneTwister(69) # to be replaced via config later
    rand_vel_noise = randn((length(Vc), n_iters))

    stoch_poles = zeros((n_iters, n_pole_vars)) # we want each soln to be a row

    @threads for i in 1:n_iters
        Vc_stochastic = Vc + vel_stds .* rand_vel_noise[:,i]
        rhs = [Vc_stochastic; zvec]

        full_soln = lhs \ rhs
        soln = full_soln[soln_idx]
        stoch_poles[i,:] = soln'
    end
    var_cov = cov(stoch_poles; dims=1)
end


function get_soln_covariance_matrix(block_matrices, lhs_fact, results, soln_idx,
                                    weighted=true)
    PvGb = block_matrices["PvGb"]
    cm = block_matrices["cm"]
    y_obs = block_matrices["Vc"]

    if weighted == true
        weights = block_matrices["weights"]
    else
        weights = ones(length(block_matrices["weights"]))
    end

    y_pred = get_pred_solution(PvGb, block_matrices["keys"],
        results["poles"]; tri_results=results["tri_slip_rates"])
    
    n, p = size(PvGb)

    e = y_pred .- y_obs
    MSE = 1 / (n - p) * e' * e
    RMSE_string = "RMSE: " * string(sqrt(MSE))
    @info RMSE_string

    if all(x -> x == 1., weights) & (length(cm) == 0)
        var_cov = make_OLS_cov(PvGb)
    elseif any(x -> x != 1., weights) & (length(cm) == 0)
        var_cov = make_WLS_cov(PvGb, weights)
    else
        var_cov = make_CWLS_cov_iter(lhs_fact, weights, y_obs, cm,
            p, soln_idx, 1000)
    end

    var = MSE * var_cov
    standard_error_vec = sqrt.(diag(var))

    SE_string = "mean standard error: " * string(
        sum(standard_error_vec) / length(standard_error_vec))
    @info SE_string

    var
end


function get_pole_uncertainties!(poles, pole_var, vel_group_keys)
    for (i, (fix, mov)) in enumerate(vel_group_keys)
        poles[(fix, mov)] = PoleCart(poles[(fix, mov)];
                                     ex=sqrt(pole_var[i * 3 - 2, i * 3 - 2]), 
                                     ey=sqrt(pole_var[i * 3 - 1, i * 3 - 1]),
                                     ez=sqrt(pole_var[i * 3, i * 3]),
                                     cxy=pole_var[i * 3 - 2, i * 3 - 1],
                                     cxz=pole_var[i * 3 - 2, i * 3],
                                     cyz=pole_var[i * 3 - 1, i * 3])
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


function predict_slip_rates(faults, poles)
    slip_rates = Oiler.Utils.get_fault_slip_rates_from_poles(
            faults, poles)

    new_faults = []

    for (i, fault) in enumerate(faults)
        push!(new_faults, 
              Fault(
                  fault.trace,
                  fault.strike,
                  fault.dip,
                  fault.dip_dir,
                  slip_rates[i][2],
                  slip_rates[i][4], # fault.extension_err,
                  slip_rates[i][1],
                  slip_rates[i][3], # fault.dextral_err,
                  slip_rates[i][5], # fault.cde,
                  fault.lsd,
                  fault.usd,
                  fault.name,
                  fault.hw,
                  fault.fw
              )
        )
    end
    new_faults
end


function get_pred_solution(PvGb, keys, poles; tri_results=Dict())
    
    soln_vec = [[poles[k].x poles[k].y poles[k].z] 
                 for k in keys]
    soln_vec = [(soln_vec...)...]
    
    if length(tri_results) > 0
        tri_soln = [[tri["dip_slip"] tri["strike_slip"]] 
                    for tri in values(tri_results)]
        tri_soln = [(tri_soln...)...]
        append!(soln_vec, tri_soln)
    end

    # multiply for pred vels
    pred_vel_vec = PvGb * soln_vec
end


function predict_model_velocities(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}},
    block_matrices, poles; tri_results=Dict())
    
    pred_vel_vec = get_pred_solution(block_matrices["PvGb"], block_matrices["keys"],
        poles; tri_results)

    # loop through vels, make new vels for each w/ predicted output
    # return in some form or fashion (pred_vel_groups?)
    pred_vels = Dict()

    counter = 1
    for pole_key in block_matrices["keys"]
        pred_vels[pole_key] = []# 

        for vel in vel_groups[pole_key]
            ve = pred_vel_vec[counter]
            counter += 1
            vn = pred_vel_vec[counter]
            counter += 1
            vu = pred_vel_vec[counter]
            counter += 1

            push!(pred_vels[pole_key], VelocityVectorSphere(vel, ve=ve, vn=vn))
        end
        convert(Array{VelocityVectorSphere}, pred_vels[pole_key])
    end
    pred_vels
end


function calc_forward_velocities(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}},
    poles::Dict{Any,Any}; faults::Array=[])

    block_inv_setup = set_up_block_inv_w_constraints(vel_groups; faults=faults,
        weighted=weighted)
    @warn "NOT FINISHED IMPLEMENTING"

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

    for iter in 1:n_iters
        Vc = build_Vc_from_vel_samples(rand_vels, vd["keys"], iter)
        rhs = make_block_inv_rhs(vd["PvGb"], Vc, zeros(p))

        kkt_soln = lhs \ rhs

        results["poles"][iter] = [PoleCart(x=kkt_soln[i * 3 - 2],
                                                y=kkt_soln[i * 3 - 1],
                                                z=kkt_soln[i * 3],
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

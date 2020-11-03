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
is zero. The latter defaults to 1e-10.
"""
function weight_from_error(error::Float64; zero_err_weight::Float64=1e-10)
    if error == 0.
        error = zero_err_weight
    end
    weight = error^-2
    # weight = error^2
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

    @info "   calculating locking effects"
    @time locking_partials = Oiler.Elastic.calc_locking_effects(faults, vel_groups)
    
    # The commented out code is a way of adding constraints to a sparse
    # matrix.  It does not seem to work, though I don't know why--
    # Sometimes the number of rows is less after the addition of more data.

    # @info "   adding to PvGb"
    # println(size(PvGb))
    # PvGb_dict = Oiler.Utils.sparse_to_dict(PvGb)
    # PvGb = [] # save some ram

    # @time for ((row_idx, col_idx), partials) in locking_partials
    #    for (i, rr) in enumerate(row_idx)
    #        for (j, cc) in enumerate(col_idx)
    #            if haskey(PvGb_dict, (rr, cc))
    #                PvGb_dict[(rr, cc)] += partials[i,j]
    #            else
    #                PvGb_dict[(rr, cc)] = partials[i,j]
    #            end
    #        end
    #    end
    # end
    # PvGb = Oiler.Utils.dict_to_sparse(PvGb_dict)
    @info "   adding to PvGb"
    PvGb = Matrix(PvGb)
    @time for (part_idx, partials) in locking_partials
        PvGb[part_idx[1], part_idx[2]] = PvGb[part_idx[1], part_idx[2]] + partials
    end
    PvGb = sparse(PvGb)
    PvGb
end


function add_tris_to_PvGb(tris, vel_groups, PvGb)
    gnss_vels = get_gnss_vels(vel_groups)
    gnss_lons = [vel["vel"].lon for vel in gnss_vels]
    gnss_lats = [vel["vel"].lat for vel in gnss_vels]
    gnss_idxs = [vel["idx"] for vel in gnss_vels]
    
    tri_effects = Oiler.Elastic.calc_tri_effects(tris, gnss_lons, gnss_lats)

    tri_eqn_matrix = zeros((size(PvGb)[1], size(tri_effects)[2]))

    for (i, vel_idx) in enumerate(gnss_idxs)
        i3 = i * 3
        gnss_row_idxs = i3 - 2:i3
        pvgb_row_idxs = vel_idx[1]

        tri_eqn_matrix[pvgb_row_idxs, :] = tri_effects[gnss_row_idxs, :]
    end

    PvGb = hcat(PvGb, tri_eqn_matrix)
end


function weight_inv_matrices(PvGb_in, Vc_in, weights)
    println(size(PvGb_in))
    N = PvGb_in' * sparse(diagm(weights))
    Vc = N * Vc_in
    PvGb = N * PvGb_in
    println(size(PvGb))
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
        lhs = [PvGb; cm .* 1e20]
        rhs = [Vc; zeros(p)]
    end

    (lhs, rhs)
end


function make_weighted_constrained_lls_matrices(PvGb, Vc, cm, weights; sparse_lhs::Bool=false)
    W = sparse(diagm(1 ./ weights))
    p, q = size(cm)
    n = length(Vc)
    if sparse_lhs
        _zeros = spzeros
    else
        _zeros = zeros
    end
    lhs = [_zeros(p, p) _zeros(p, n) cm;
           _zeros(n, p) W           PvGb;
           cm'           PvGb'       _zeros(q, q)]

    rhs = [zeros(p); Vc; zeros(q)]


    if sparse_lhs
        return sparse(lhs), rhs
    else
        return lhs, rhs
    end
end

function set_up_block_inv_w_constraints(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}};
    faults::Array=[], tris::Array=[], weighted::Bool=true, regularize::Bool=false,
    l2_lambda::Float64=100.0, sparse_lhs::Bool=false)

    @info " making block inversion matrices"
    vd = make_block_inversion_matrices_from_vels(vel_groups)
    PvGb_size = size(vd["PvGb"])
    @info "Raw PvGb size: $PvGb_size"

    @info " done making block inversion matrices"
    @info " finding vel cycles"
    @time cycles = find_vel_cycles(vd["keys"])
    @info " done finding vel cycles"

    if length(faults) > 0
        @info " doing locking"
        PvGb = add_fault_locking_to_PvGb(faults, vel_groups, vd["PvGb"])
        @info " done doing locking"
    else
        PvGb = vd["PvGb"]
    end

    if length(tris) > 0
        @info " doing tris"
        PvGb = add_tris_to_PvGb(tris, vel_groups, PvGb)
    end

    Vc = vd["Vc"]
    weights = vd["weights"]

    if cycles == Dict{Any,Any}()
        # lhs = PvGb
        # rhs = Vc
        cm = []
    else
        cm = build_constraint_matrices(cycles, vd["keys"])
        if length(tris) > 0
            cm = hcat(cm, zeros(size(cm)[1], length(tris) * 2))
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
                weights; sparse_lhs=sparse_lhs)
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

    Dict("lhs" => lhs, "rhs" => rhs, "keys" => vd["keys"], 
         "n_constraints" => size(cm, 1))
end


function solve_block_invs_from_vel_groups(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}};
    faults::Array=[], tris=[], weighted::Bool=true, regularize::Bool=false,
    l2_lambda::Float64=100.0, check_closures::Bool=true, sparse_lhs::Bool=false)

    @info "setting up matrices"
    @time block_inv_setup = set_up_block_inv_w_constraints(vel_groups; 
        faults=faults, 
        tris=tris,
        weighted=weighted, 
        regularize=regularize, l2_lambda=l2_lambda,
        sparse_lhs=sparse_lhs)
    
    lhs = block_inv_setup["lhs"]
    lhs_size = size(lhs)
    rhs_size = size(block_inv_setup["rhs"])
    @info "LHS block size $lhs_size"
    @info "RHS size $rhs_size"

    @info "solving"
    @time kkt_soln = lhs \ block_inv_setup["rhs"]

    nk = length(kkt_soln)
    np = length(block_inv_setup["keys"])
    nt = length(tris) * 2

    if weighted == true
        soln_idx = ((nk - nt) - np * 3) + 1:(nk - nt)
    else
        soln_idx = 1:np * 3
    end
    soln = kkt_soln[soln_idx]
    if nt > 0
        tri_results = Dict()
        # println(nk)
        # println(nt)
        tri_soln_idx = last(soln_idx) + 1:last(soln_idx) + nt
        # println(tri_soln_idx)
        tri_soln = kkt_soln[tri_soln_idx]
        # println(tri_soln)

        for (i, tri) in enumerate(tris)
            ds_ind = i * 2 - 1
            ss_ind = i * 2
            tri_results[tri.name] = Dict()
            tri_results[tri.name]["dip_slip"] = tri_soln[ds_ind]
            tri_results[tri.name]["strike_slip"] = tri_soln[ss_ind]
        end

    end


    poles = Dict()
    for (i, (fix, mov)) in enumerate(block_inv_setup["keys"])
        poles[(fix, mov)] = PoleCart(x=soln[i * 3 - 2], 
                                     y=soln[i * 3 - 1],
                                     z=soln[i * 3],
                                     fix=fix, mov=mov)
    end

    if check_closures == true
        closures = Oiler.Utils.check_vel_closures(poles)
    end

    if nt > 0
        return poles, tri_results
    else
        return poles
    end
end


function calc_forward_velocities(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}},
    poles::Dict{Any,Any}; faults::Array=[], weighted::Bool=true)

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

module Strain

export add_strain_to_PvGb, estimate_block_strain_rates,
    get_principal_orientation, get_principal_stresses,
    get_block_strain_results, get_block_strain_uncertainties!

using ..Oiler

using Logging
using SparseArrays
using LinearAlgebra
using DataFrames

const M_TO_MM = 1e3
const NSTRAIN_TO_STRAIN = 1e-9


function normalize_strain_block_id(val)
    if ismissing(val) || isnothing(val)
        return ""
    end
    string(strip(string(val)))
end


function normalize_strain_block_ids(block_ids)
    out = String[]
    for block_id in block_ids
        block_id = normalize_strain_block_id(block_id)
        if !isempty(block_id) && !(block_id in out)
            push!(out, block_id)
        end
    end
    out
end


function get_block_reference_points(block_df; epsg=102016)
    ref_points = Dict{String,Tuple{Float64,Float64}}()
    block_ids = string.(block_df[!, :fid])
    colnames = Symbol.(names(block_df))

    if (:lon in colnames) && (:lat in colnames)
        lons = Float64.(block_df[!, :lon])
        lats = Float64.(block_df[!, :lat])
    elseif :geometry in colnames
        centroids = Oiler.Utils.get_block_centroids(block_df; epsg=epsg)
        lons = [centroid.coords[1] for centroid in centroids]
        lats = [centroid.coords[2] for centroid in centroids]
    else
        err_msg = "block_df must contain either lon/lat columns or geometry for strain calculations"
        throw(ArgumentError(err_msg))
    end

    for (i, block_id) in enumerate(block_ids)
        ref_points[block_id] = (lons[i], lats[i])
    end
    ref_points
end


function get_local_xy_m(lon::Float64, lat::Float64, ref_lon::Float64, ref_lat::Float64)
    x, y = Oiler.Geom.azimuthal_equidistant_proj([lon], [lat], ref_lon, ref_lat)
    x[1], y[1]
end


function make_strain_rate_equations(lon::Float64, lat::Float64, ref_lon::Float64,
    ref_lat::Float64; include_translation::Bool=false)

    x_m, y_m = get_local_xy_m(lon, lat, ref_lon, ref_lat)
    x_coeff = x_m * M_TO_MM * NSTRAIN_TO_STRAIN
    y_coeff = y_m * M_TO_MM * NSTRAIN_TO_STRAIN

    if include_translation
        lhs = [1.0 0.0 x_coeff y_coeff 0.0;
               0.0 1.0 0.0 x_coeff y_coeff]
    else
        lhs = [x_coeff y_coeff 0.0;
               0.0 x_coeff y_coeff]
    end
    lhs
end


function count_gnss_data_per_block(vel_groups)
    counts = Dict{String,Int}()
    for (_, vels) in vel_groups
        for vel in vels
            if vel.vel_type == "GNSS"
                mov = normalize_strain_block_id(vel.mov)
                if !isempty(mov)
                    counts[mov] = get(counts, mov, 0) + 1
                end
            end
        end
    end
    counts
end


function get_active_strain_blocks(vel_groups, block_df, strain_blocks; min_data::Int=5,
    epsg=102016)
    ref_points = get_block_reference_points(block_df; epsg=epsg)
    if isempty(strain_blocks)
        requested_blocks = string.(block_df[!, :fid])
    else
        requested_blocks = normalize_strain_block_ids(strain_blocks)
    end
    data_counts = count_gnss_data_per_block(vel_groups)

    active_blocks = String[]
    for block_id in requested_blocks
        if !haskey(ref_points, block_id)
            @warn "strain block $block_id not found in block_df, skipping"
            continue
        end
        n_obs = get(data_counts, block_id, 0)
        if n_obs < min_data
            @info "skipping strain block $block_id: $n_obs GNSS-type observations < min_data=$min_data"
            continue
        end
        push!(active_blocks, block_id)
    end

    active_blocks, ref_points, data_counts
end


function build_block_strain_matrix(vel_groups, vel_group_keys, block_df, strain_blocks;
    min_data::Int=5, epsg=102016)
    active_blocks, ref_points, data_counts = get_active_strain_blocks(
        vel_groups, block_df, strain_blocks; min_data=min_data, epsg=epsg)

    n_rows = sum(3 * length(vel_groups[key]) for key in vel_group_keys)
    n_params = 3 * length(active_blocks)
    if n_params == 0
        return spzeros(n_rows, 0), active_blocks, ref_points, data_counts
    end

    block_col_start = Dict(block_id => 3 * (i - 1) + 1 for (i, block_id) in enumerate(active_blocks))
    row_idxs = Int[]
    col_idxs = Int[]
    vals = Float64[]

    row_start = 1
    for key in vel_group_keys
        for vel in vel_groups[key]
            if vel.vel_type == "GNSS" && haskey(block_col_start, vel.mov)
                ref_lon, ref_lat = ref_points[vel.mov]
                strain_rows = make_strain_rate_equations(vel.lon, vel.lat, ref_lon, ref_lat)
                col_start = block_col_start[vel.mov]
                for local_row in 1:2
                    for local_col in 1:3
                        val = strain_rows[local_row, local_col]
                        if val != 0.0
                            push!(row_idxs, row_start + local_row - 1)
                            push!(col_idxs, col_start + local_col - 1)
                            push!(vals, val)
                        end
                    end
                end
            end
            row_start += 3
        end
    end

    strain_mat = sparse(row_idxs, col_idxs, vals, n_rows, n_params)
    strain_mat, active_blocks, ref_points, data_counts
end


function make_block_strain_prior_matrices(strain_block_ids; strain_rate_err::Float64=5.0)
    n_params = 3 * length(strain_block_ids)
    if n_params == 0
        return zeros(0, 0), zeros(0), spzeros(0, 0)
    end
    if !(strain_rate_err > 0.0) || !isfinite(strain_rate_err)
        throw(ArgumentError("strain_rate_err must be finite and > 0"))
    end

    lhs = Matrix{Float64}(I, n_params, n_params)
    rhs = zeros(n_params)
    weights = spdiagm(0 => fill(strain_rate_err^2, n_params))
    lhs, rhs, weights
end


function add_strain_to_PvGb(vel_groups, vd, block_df, strain_blocks;
    regularize::Bool=true, strain_rate_err::Float64=5.0,
    min_data::Int=5, epsg=102016)

    if isnothing(block_df)
        throw(ArgumentError("block_df is required when solving for intra-block strain"))
    end

    PvGb = vd["PvGb"]
    PvGb_n_cols = size(PvGb, 2)
    strain_mat, strain_block_ids, ref_points, data_counts = build_block_strain_matrix(
        vel_groups, vd["keys"], block_df, strain_blocks; min_data=min_data, epsg=epsg)

    vd["strain_block_ids"] = strain_block_ids
    vd["strain_center_lons"] = [ref_points[block_id][1] for block_id in strain_block_ids]
    vd["strain_center_lats"] = [ref_points[block_id][2] for block_id in strain_block_ids]
    vd["strain_data_counts"] = [get(data_counts, block_id, 0) for block_id in strain_block_ids]

    if size(strain_mat, 2) == 0
        return vd
    end

    replace!(strain_mat.nzval, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
    PvGb = hcat(PvGb, strain_mat)

    if regularize
        strain_reg_matrix, strain_reg_vels, strain_reg_weights = make_block_strain_prior_matrices(
            strain_block_ids; strain_rate_err=strain_rate_err)
        strain_reg_block = hcat(spzeros(size(strain_reg_matrix, 1), PvGb_n_cols),
            sparse(strain_reg_matrix))
        PvGb = vcat(PvGb, strain_reg_block)
        vd["Vc"] = vcat(vd["Vc"], strain_reg_vels)
        vd["extra_weights"] = Oiler.Utils.diagonalize_matrices(
            [get(vd, "extra_weights", spzeros(0, 0)), strain_reg_weights])
        vd["strain_weights"] = strain_reg_weights
    end

    dropzeros!(PvGb)
    vd["PvGb"] = PvGb
    vd
end


function get_block_strain_results(strain_soln, strain_block_ids)
    results = Dict{String,Dict{String,Float64}}()
    for (i, block_id) in enumerate(strain_block_ids)
        row_start = 3 * (i - 1) + 1
        row_range = row_start:(row_start + 2)
        results[block_id] = Dict(
            "ee" => strain_soln[row_range[1]],
            "en" => strain_soln[row_range[2]],
            "nn" => strain_soln[row_range[3]],
        )
    end
    results
end


function get_block_strain_uncertainties!(soln_cov, strain_results, strain_soln_idx,
    strain_block_ids)
    if length(strain_block_ids) == 0
        return
    end

    strain_cov = soln_cov[strain_soln_idx, strain_soln_idx]
    for (i, block_id) in enumerate(strain_block_ids)
        row_start = 3 * (i - 1) + 1
        row_range = row_start:(row_start + 2)
        block_cov = Matrix{Float64}(strain_cov[row_range, row_range])
        strain_results[block_id]["ee_err"] = sqrt(max(block_cov[1, 1], 0.0))
        strain_results[block_id]["en_err"] = sqrt(max(block_cov[2, 2], 0.0))
        strain_results[block_id]["nn_err"] = sqrt(max(block_cov[3, 3], 0.0))
        strain_results[block_id]["ee_en_cov"] = block_cov[1, 2]
        strain_results[block_id]["ee_nn_cov"] = block_cov[1, 3]
        strain_results[block_id]["en_nn_cov"] = block_cov[2, 3]
    end
end


function estimate_block_strain_rates(block_id, vel_df, block_df; min_stations::Int=5,
    ve=:ve, vn=:vn, lon=:lon, lat=:lat, epsg=102016,
    fit_translation::Bool=true)

    block_id = normalize_strain_block_id(block_id)
    vels = vel_df[(string.(vel_df.mov) .== block_id), :]

    if size(vels, 1) < min_stations
        return (ee=0.0, en=0.0, nn=0.0, ve0=0.0, vn0=0.0)
    end

    ref_points = get_block_reference_points(block_df; epsg=epsg)
    if !haskey(ref_points, block_id)
        throw(ArgumentError("block $block_id not found in block_df"))
    end
    ref_lon, ref_lat = ref_points[block_id]

    lhs_rows = Matrix{Float64}[]
    rhs = Float64[]
    for row in eachrow(vels)
        lhs = make_strain_rate_equations(Float64(row[lon]), Float64(row[lat]),
            ref_lon, ref_lat; include_translation=fit_translation)
        push!(lhs_rows, lhs)
        push!(rhs, Float64(row[ve]))
        push!(rhs, Float64(row[vn]))
    end

    lhs = reduce(vcat, lhs_rows)
    soln = lhs \ rhs

    if fit_translation
        return (ee=soln[3], en=soln[4], nn=soln[5], ve0=soln[1], vn0=soln[2])
    else
        return (ee=soln[1], en=soln[2], nn=soln[3], ve0=0.0, vn0=0.0)
    end
end


function get_principal_orientation(E)
    E = Symmetric(E)
    if iszero(norm(E))
        return 0.0
    end
    0.5 * atan(2 * E[1, 2], E[1, 1] - E[2, 2])
end


function get_principal_stresses(E)
    sort(eigvals(Symmetric(E)), rev=true)
end

end

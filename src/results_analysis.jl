module ResultsAnalysis

using Setfield
using StatsBase# , DataAPI
import Base.Threads.@threads
using DataFrames, DataFramesMeta

using ..Oiler

function check_missing_rate(val)
    if ismissing(val) | isnothing(val)
        return false
    elseif val == ""
        return false
    elseif val == 0.0
        return false
    elseif ismissing(val) | isnothing(val)
        return false
    else
        return true
    end
end


function get_geol_slip_rate_indices(geol_slip_rate_df)
    include_idx = map(!, iszero.(geol_slip_rate_df[!, :include]))

    dex_has_rate_idx = map(check_missing_rate,
        geol_slip_rate_df[!, :dextral_rate])
    ext_has_rate_idx = map(check_missing_rate,
        geol_slip_rate_df[!, :extension_rate])

    include_idx, dex_has_rate_idx, ext_has_rate_idx
end


function get_geol_obs_slip_rates(geol_slip_rate_df)

    include_idx, dex_has_rate_idx, ext_has_rate_idx = get_geol_slip_rate_indices(
        geol_slip_rate_df
    )

    dex_include = include_idx .* dex_has_rate_idx
    ext_include = include_idx .* ext_has_rate_idx

    dex_geol_obs = geol_slip_rate_df[!, :dextral_rate][dex_include]
    dex_geol_err = geol_slip_rate_df[!, :dextral_err][dex_include]

    ext_geol_obs = geol_slip_rate_df[!, :extension_rate][ext_include]
    ext_geol_err = geol_slip_rate_df[!, :extension_err][ext_include]

    Dict("dex_geol_obs" => dex_geol_obs,
        "dex_geol_err" => dex_geol_err,
        "ext_geol_obs" => ext_geol_obs,
        "ext_geol_err" => ext_geol_err)
end


function get_geol_pred_slip_rates(geol_slip_rate_vels, fault_df, results;
    usd=:usd, lsd=:lsd)
    pred_geol_slip_rates = []
    for (i, rate) in enumerate(geol_slip_rate_vels)
        fault_idx = rate.name
        fault_row = @subset(fault_df, :fid .== fault_idx)[1, :]
        fault = Oiler.IO.row_to_fault(fault_row; lsd=lsd,
            usd=usd)

        if haskey(results["poles"], (rate.fix, rate.mov))
            pred_rate = Oiler.Faults.get_fault_slip_rate_from_pole(fault,
                results["poles"][(rate.fix, rate.mov)];
                lon=rate.lon, lat=rate.lat)
        else
            pred_rate = Oiler.Faults.get_fault_slip_rate_from_pole(fault,
                results["poles"][(rate.mov, rate.fix)];
                lon=rate.lon, lat=rate.lat)
        end

        push!(pred_geol_slip_rates, pred_rate)
    end
    pred_geol_slip_rates
end


function make_geol_slip_rate_obs_res_table(geol_slip_rate_vels, geol_slip_rate_df, fault_df, results;
    usd=:usd, lsd=:lsd)

    pr = get_geol_pred_slip_rates(geol_slip_rate_vels, fault_df, results; usd=usd, lsd=lsd)

    v_rl_pred = collect(p[1] for p in pr)
    e_rl_pred = collect(p[3] for p in pr)
    v_ex_pred = collect(p[2] for p in pr)
    e_ex_pred = collect(p[4] for p in pr)
    cde_pred  = collect(p[5] for p in pr)

    new_rate_df = copy(geol_slip_rate_df)

    new_rate_df[!, :dextral_rate_pred] = v_rl_pred
    new_rate_df[!, :dextral_err_pred] = e_rl_pred
    new_rate_df[!, :extension_rate_pred] = v_ex_pred
    new_rate_df[!, :extension_err_pred] = e_ex_pred
    new_rate_df[!, :cde_pred] = cde_pred

    new_rate_df[!,:lon], new_rate_df[!, :lat] = Oiler.Geom.get_coord_vecs(geol_slip_rate_df[!, :geometry])
    
    new_rate_df = new_rate_df[!,
        [:dextral_rate, :dextral_err, :dextral_rate_pred, :dextral_err_pred, 
        :extension_rate, :extension_err, :extension_rate_pred, :extension_err_pred, 
        :cde_pred, :lon, :lat, :source
        ]]
end


function filter_geol_pred_rates(geol_pred_slip_rates, include_idx,
    dex_has_rate_idx, ext_has_rate_idx)

    dex_inc = dex_has_rate_idx[include_idx]
    ext_inc = ext_has_rate_idx[include_idx]

    dex_geol_pred = [p[1] for p in geol_pred_slip_rates[dex_inc]]
    dex_geol_err = [p[3] for p in geol_pred_slip_rates[dex_inc]]

    ext_geol_pred = [p[2] for p in geol_pred_slip_rates[ext_inc]]
    ext_geol_err = [p[4] for p in geol_pred_slip_rates[ext_inc]]


    Dict("dex_geol_pred" => dex_geol_pred,
        "dex_geol_pred_err" => dex_geol_err,
        "ext_geol_pred" => ext_geol_pred,
        "ext_geol_pred_err" => ext_geol_err)
end


function get_obs_pred_geol_rate_vecs(; geol_slip_rate_df,
    geol_slip_rate_vels,
    fault_df,
    results,
    usd=:usd, lsd=:lsd)

    include_idx, dex_has_rate_idx, ext_has_rate_idx = get_geol_slip_rate_indices(
        geol_slip_rate_df
    )

    pred_rates = get_geol_pred_slip_rates(geol_slip_rate_vels, fault_df,
        results; usd=usd, lsd=lsd)

    pred_rates = filter_geol_pred_rates(pred_rates, include_idx,
        dex_has_rate_idx, ext_has_rate_idx)

    obs_rates = get_geol_obs_slip_rates(geol_slip_rate_df)

    obs_rates, pred_rates
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


function calc_RMSE_from_G(block_matrices, results)
    y_pred = get_pred_solution(block_matrices["PvGb"], block_matrices["keys"],
        results["poles"];
        tri_results=results["tri_slip_rates"])

    y_obs = block_matrices["Vc"]
    n, p = size(block_matrices["PvGb"])

    RMSE = Oiler.Stats.RMSE(y_obs, y_pred, n, p)
end

function predict_model_velocities(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}},
    block_matrices, poles; tri_results=Dict())

    # this is for propagating uncertainties to GNSS vels, unfinished
    # but I don't want to leave it in a branch.
    #gnss_df = Oiler.Utils.make_gnss_df_from_vel_groups(vel_groups)

    #if haskey(block_matrices, "stoch_poles")
    #    n_iters = size(block_matrices["stoch_poles"], 1)

    #    #pred_vels_stoch = Dict()
    #    gnss_stoch_ve = zeros(size(gnss_df, 1), n_iters)
    #    gnss_stoch_vn = zeros(size(gnss_df, 1), n_iters)

    #    @warn "doing the thing"
    #    @time @threads for i in 1:n_iters
    #        poles_iter = Oiler.Solver.get_poles_from_soln(block_matrices["keys"],
    #            block_matrices["stoch_poles"][i, :])

    #        # I don't think I'm propagating tri uncertainties through this...
    #        pred_vels_stoch_i = predict_model_velocities_one_pole_set(vel_groups,
    #            block_matrices, poles_iter; tri_results=tri_results)

    #        pred_vels = [v["vel"] for v in Oiler.Utils.get_gnss_vels(
    #            pred_vels_stoch_i)]
    #        pve = [v.ve for v in pred_vels]
    #        pvn = [v.vn for v in pred_vels]

    #        gnss_stoch_ve[:, i] = pve
    #        gnss_stoch_vn[:, i] = pvn
    #    end
    #end
    predict_model_velocities_one_pole_set(vel_groups,
        block_matrices, poles; tri_results=tri_results)
end


function predict_model_velocities_one_pole_set(vel_groups::Dict{Tuple{String,String},Array{VelocityVectorSphere,1}},
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
                fault.fw,
                fault.fid
            )
        )
    end
    new_faults
end


function _get_obs_pred_vecs(; results,
    vel_groups,
    geol_slip_rate_df,
    geol_slip_rate_vels,
    fault_df,
    usd=:usd, lsd=:lsd
)

    geol_data_obs, geol_data_pred = get_obs_pred_geol_rate_vecs(
        geol_slip_rate_df=geol_slip_rate_df,
        geol_slip_rate_vels=geol_slip_rate_vels,
        fault_df=fault_df,
        results=results,
        usd=usd, lsd=lsd)

    gnss_results = get_gnss_results(results, vel_groups)

    obs_vec = vcat(gnss_results.obs_ve, gnss_results.obs_vn,
        geol_data_obs["dex_geol_obs"], geol_data_obs["ext_geol_obs"])
    obs_err_vec = vcat(gnss_results.obs_ee, gnss_results.obs_en,
        geol_data_obs["dex_geol_err"], geol_data_obs["ext_geol_err"])
    pred_vec = vcat(gnss_results.pred_ve, gnss_results.pred_vn,
        geol_data_pred["dex_geol_pred"], geol_data_pred["ext_geol_pred"])
    #pred_err_vec = vcat(gnss_results.pred_ee, gnss_results.pred_en, 
    #               geol_data_pred["dex_geol_pred_err"], geol_data_pred["ext_geol_pred_err"])

    out = Dict("geol_obs"=>geol_data_obs, "geol_pred"=>geol_data_pred, "gnss"=>gnss_results)
end


function compare_data_results(; results,
    vel_groups,
    geol_slip_rate_df,
    geol_slip_rate_vels,
    fault_df,
    usd=:usd, lsd=:lsd
)

    geol_data_obs, geol_data_pred = get_obs_pred_geol_rate_vecs(
        geol_slip_rate_df=geol_slip_rate_df,
        geol_slip_rate_vels=geol_slip_rate_vels,
        fault_df=fault_df,
        results=results,
        usd=usd, lsd=lsd)

    geol_rates_table = make_geol_slip_rate_obs_res_table(geol_slip_rate_vels, 
        geol_slip_rate_df, fault_df, results, usd=usd, lsd=lsd)
    gnss_results = get_gnss_results(results, vel_groups)

    obs_vec = vcat(gnss_results.obs_ve, gnss_results.obs_vn,
        geol_data_obs["dex_geol_obs"], geol_data_obs["ext_geol_obs"])
    obs_err_vec = vcat(gnss_results.obs_ee, gnss_results.obs_en,
        geol_data_obs["dex_geol_err"], geol_data_obs["ext_geol_err"])
    pred_vec = vcat(gnss_results.pred_ve, gnss_results.pred_vn,
        geol_data_pred["dex_geol_pred"], geol_data_pred["ext_geol_pred"])


    # pred_err_vec = vcat(gnss_results.pred_ee, gnss_results.pred_en, 
    #               geol_data_pred["dex_geol_pred_err"], geol_data_pred["ext_geol_pred_err"])

    resid_vec = obs_vec .- pred_vec

    # results["stats_info"]["summary_stats"]

    gnss_misfits = gnss_chi_sq_misfits(gnss_results)
    geol_misfits = geol_chi_sq_misfits(obs_vec, pred_vec, obs_err_vec)

    misfits = vcat(gnss_misfits, geol_misfits)

    results["stats_info"]["reduced_chi_sq"] = Oiler.Stats.reduced_chi_sq_mccaffrey(
        misfits, results["stats_info"]["n_params"]
    )

    resid_summary = StatsBase.summarystats(abs.(resid_vec))

    results["data"] = Dict(
        "geol_data_obs" => geol_data_obs,
        "geol_data_pred" => geol_data_pred,
        "gnss_results" => gnss_results,
        "geol_rates_table" => geol_rates_table,
    )

    results["stats_info"]["resid_summary"] = Dict(
        "mean" => resid_summary.mean,
        "min" => resid_summary.min,
        "q25" => resid_summary.q25,
        "median" => resid_summary.median,
        "q75" => resid_summary.q75,
        "max" => resid_summary.max,
    )

end


function get_gnss_results(results, vel_groups)
    get_gnss_results_per_iter(results["predicted_vels"], vel_groups)
end


function gnss_chi_sq_misfits(gnss_results)
    misfits = map(x->Oiler.Stats.vector_misfit(
        x.obs_ve, x.obs_vn, x.pred_ve, x.pred_vn,
        x.obs_ee, x.obs_en, x.cen),
        eachrow(gnss_results))
end
    

function geol_chi_sq_misfits(geol_obs, geol_pred, geol_err)
    return (geol_obs .- geol_pred).^2 ./ geol_err .^ 2
end


function get_gnss_results_per_iter(predicted_vels, vel_groups)

    gnss_df = Oiler.Utils.make_gnss_df_from_vel_groups(vel_groups)
    rename!(gnss_df, :ve => :obs_ve)
    rename!(gnss_df, :vn => :obs_vn)
    rename!(gnss_df, :ee => :obs_ee)
    rename!(gnss_df, :en => :obs_en)


    pred_vels = [v["vel"] for v in Oiler.Utils.get_gnss_vels(
        predicted_vels)]
    pve = [v.ve for v in pred_vels]
    pvn = [v.vn for v in pred_vels]

    rve, rvn = gnss_df.obs_ve - pve, gnss_df.obs_vn - pvn

    gnss_df.pred_ve = round.(pve, digits=3)
    gnss_df.pred_vn = round.(pvn, digits=3)
    # these are not correctly estimated
    # gnss_df.pred_ee =  round.([v.ee for v in pred_vels], digits=3)
    # gnss_df.pred_en =  round.([v.en for v in pred_vels], digits=3)
    gnss_df.cen = round.([v.cen for v in pred_vels], digits=3)
    gnss_df.resid_e = round.(rve, digits=3)
    gnss_df.resid_n = round.(rvn, digits=3)
    # gnss_df.resid_ee = round.(sqrt.(gnss_df.pred_ee.^2 + gnss_df.obs_ee.^2),
    #                               digits=3)
    # gnss_df.resid_en = round.(sqrt.(gnss_df.pred_en.^2 + gnss_df.obs_en.^2),
    #                          digits=3)
    gnss_df
end



function get_tri_locking_rate(tri, poles; set=false)
    # if haskey(poles, (tri.fw, tri.hw))
    #    pole = poles[(tri.fw, tri.hw)]
    # elseif haskey(poles, (tri.hw, tri.fw))
    #    pole = -poles[(tri.hw, tri.fw)]
    # else
    # try
    pole = Oiler.Utils.get_path_euler_pole(poles, tri.hw, tri.fw)
    # catch
    #    tri_name = tri.name
    #    tri_pole = (tri.fw, tri.hw)
    #    warn_msg = "$tri_name pole $tri_pole not in pole results"
    #    @warn warn_msg
    #    return missing
    # end

    ds, de, ss, se, cds = Oiler.Tris.get_tri_rate_from_pole(tri, pole)
    dip_locking_frac = ds / tri.dip_slip_rate
    strike_locking_frac = ss / tri.strike_slip_rate

    if set
        tri = @set tri.dip_locking_frac = dip_locking_frac
        tri = @set tri.strike_locking_frac = strike_locking_frac
        return tri
    else
        return dip_locking_frac, strike_locking_frac
    end
end


function get_block_centroid_vels(results, block_df; fix)
    centroid_df = Oiler.Utils.get_block_centroid_velocities(block_df,
        results["poles"]; fix=fix)
    results["block_centroids"] = centroid_df
    return
end


function calculate_resid_block_strain_rates(results; min_stations=5.)
    if !(haskey(results, "block_centroids"))
        error("Must calculate block centroids first")
    end

    blocks = results["block_centroids"]
    vels = results["data"]["gnss_results"]

    n_blocks = size(blocks, 1)

    strain_ee = zeros(n_blocks)
    strain_en = zeros(n_blocks)
    strain_nn = zeros(n_blocks)
    mc_az = zeros(n_blocks)
    mc_stress = zeros(n_blocks)
    lc_stress = zeros(n_blocks)

    for (i, block) in enumerate(eachrow(blocks))
        block_strain = Oiler.Strain.get_block_strain_rates(block.fid, vels, blocks; 
                                              min_stations=min_stations,
            ve=:resid_e, vn=:resid_n)
    
        strain_ee[i] = block_strain.ee
        strain_en[i] = block_strain.en
        strain_nn[i] = block_strain.nn

        E = [block_strain.ee block_strain.en; block_strain.en block_strain.nn]
        mc_az[i] = Oiler.Geom.angle_to_az(
            Oiler.Strain.get_principal_orientation(E)
        )

        princ_stresses = Oiler.Strain.get_principal_stresses(E)
        mc_stress[i] = princ_stresses[1]
        lc_stress[i] = princ_stresses[2]
    end

    blocks[!,"mc_stress"] = mc_stress
    blocks[!,"lc_stress"] = lc_stress
    blocks[!,"mc_azimuth"] = mc_az
    return
end





end # module

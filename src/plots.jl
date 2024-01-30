module Plots

using PyPlot
using DataFrames, DataFramesMeta
import Statistics: quantile

using Oiler


function plot_tri_prior_post(tris, results)

    fig = figure()

    prior_rates = [sqrt(tri.dip_slip_rate^2 + tri.strike_slip_rate^2) for tri in tris]
    post_rates = [get_tri_total_rate(tri, results) for tri in tris]

    plot(prior_rates, post_rates, ".")

    data_min = min(minimum(prior_rates), minimum(post_rates))
    data_max = max(maximum(prior_rates), maximum(post_rates))

    plot([data_min, data_max], [data_min, data_max], "C1--", lw=0.5)

    axis("equal")

    xlabel("prior rate")
    ylabel("posterior rate")

    title("tri results")

    return fig
end


function plot_results_map(results, vel_groups, faults, tris=[])
    vel_df = Oiler.ResultsAnalysis.get_gnss_results(results, vel_groups)

    fig, ax = subplots(figsize=(14, 14))

    for fault in faults
        plot(fault.trace[:, 1], fault.trace[:, 2], "k-", lw=0.5)
    end

    if length(tris) > 0
        cm = get_cmap(:viridis)

        tri_rates = [get_tri_total_rate(tri, results) for tri in tris]
        tri_rate_min = minimum(tri_rates)
        tri_rate_max = maximum(tri_rates)

        #clim(tri_rate_min, tri_rate_max)
        norm = PyPlot.matplotlib.colors.Normalize(tri_rate_min, tri_rate_max)
        mappable = PyPlot.matplotlib.cm.ScalarMappable(norm=norm, cmap=cm)

        for tri in tris
            plot_tri(tri, results; vmin=tri_rate_min, vmax=tri_rate_max, cm=cm)
        end
        colorbar(mappable, ax=ax)
    end

    x_limits = xlim()
    y_limits = ylim()

    quiver(vel_df.lon, vel_df.lat,
        vel_df.obs_ve, vel_df.obs_vn,
        color="b", scale=300, label="obs GNSS")
    quiver(vel_df.lon, vel_df.lat,
        vel_df.pred_ve, vel_df.pred_vn,
        color="r", scale=300, label="pred GNSS")

    axis("equal")
    xlim(x_limits)
    ylim(y_limits)
    legend()

    return fig
end

function get_tri_total_rate(tri, results)
    ds = results["tri_slip_rates"][tri.name]["dip_slip"]
    ss = results["tri_slip_rates"][tri.name]["strike_slip"]
    total_rate = sqrt(ds^2 + ss^2)
end


function plot_tri(tri, results; vmin=:tri_rate_min, vmax=:tri_rate_max, cm=:cm)
    lons = [tri.p1[1], tri.p2[1], tri.p3[1], tri.p1[1]]
    lats = [tri.p1[2], tri.p2[2], tri.p3[2], tri.p1[2]]
    total_rate = get_tri_total_rate(tri, results)
    rate_frac = (total_rate - vmin) / (vmax - vmin)
    color = cm(rate_frac)
    fill(lons, lats, color=color, alpha=0.25, zorder=0)
end


function plot_slip_rate_fig(geol_slip_rate_df, geol_slip_rate_vels,
    fault_df, results; usd=:usd, lsd=:lsd)
    pred_geol_slip_rates = []
    for (i, rate) in enumerate(geol_slip_rate_vels)
        fault_fid_type = typeof(fault_df[1, :fid])

        if fault_fid_type <: AbstractString
            fault_idx = rate.name
        else
            fault_idx = parse(fault_fid_type, rate.name)
        end
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


    inc = map(!, iszero.(geol_slip_rate_df[!, :include]))

    function check_missing(val)
        if isnothing(val) | ismissing(val) | (typeof(val) == Missing)
            return false
        elseif val == ""
            return false
        elseif val == 0.0
            return false
        else
            return true
        end
    end

    dex_inc = [check_missing(d) for d in geol_slip_rate_df[!, :dextral_rate]]
    ext_inc = [check_missing(d) for d in geol_slip_rate_df[!, :extension_rate]]

    dex_geol_obs = geol_slip_rate_df[!, :dextral_rate][dex_inc.*inc]
    # dex_geol_err = parse.(Float64, geol_slip_rate_df[!,:dextral_err][dex_inc .* inc])
    dex_geol_err = geol_slip_rate_df[!, :dextral_err][dex_inc.*inc]
    dex_geol_pred = [p[1] for p in pred_geol_slip_rates[dex_inc[inc]]]
    dex_geol_pred_err = [p[3] for p in pred_geol_slip_rates[dex_inc[inc]]]

    ext_geol_obs = geol_slip_rate_df[!, :extension_rate][ext_inc.*inc]
    ext_geol_err = geol_slip_rate_df[!, :extension_err][ext_inc.*inc]
    ext_geol_pred = [p[2] for p in pred_geol_slip_rates[ext_inc[inc]]]
    ext_geol_pred_err = [p[4] for p in pred_geol_slip_rates[ext_inc[inc]]]

    fig = figure(figsize=(5, 9))
    # suptitle("Observed vs. Modeled Quaternary Slip Rates")

    subplot(2, 1, 1)
    data_max = maximum([maximum(dex_geol_obs), maximum(dex_geol_pred)])
    data_min = minimum([minimum(dex_geol_obs), minimum(dex_geol_pred)])

    plot([data_min, data_max], [data_min, data_max], "C1--", lw=0.5)

    axis("equal")

    dex_geol_obs = convert(Array{Float64}, dex_geol_obs)

    errorbar(dex_geol_obs, dex_geol_pred, xerr=dex_geol_err, yerr=dex_geol_pred_err,
        fmt=",", elinewidth=0.3)

    autoscale(false)

    fill_between([data_min, 0.0, -data_min],
        [data_min, 0.0, data_min],
        [data_min, data_min, data_min],
        color="grey",
        lw=0.0,
        alpha=0.1)

    fill_between([data_min, 0.0, -data_min],
        [-data_min, -data_min, -data_min],
        [-data_min, 0.0, -data_min],
        color="grey",
        lw=0.0,
        alpha=0.1)

    xlabel("observed slip rate (mm a⁻¹)")
    ylabel("modeled slip rate (mm a⁻¹)")
    title("dextral")
    subplot(2, 1, 2)

    data_max = maximum([maximum(ext_geol_obs), maximum(ext_geol_pred)])
    data_min = minimum([minimum(ext_geol_obs), minimum(ext_geol_pred)])
    plot([data_min, data_max], [data_min, data_max], "C1--", lw=0.5)

    axis("equal")
    errorbar(ext_geol_obs, ext_geol_pred, xerr=ext_geol_err, yerr=ext_geol_pred_err,
        fmt=",", elinewidth=0.3)

    autoscale(false)

    fill_between([data_min, 0.0, -data_min],
        [data_min, 0.0, data_min],
        [data_min, data_min, data_min],
        color="grey",
        lw=0.0,
        alpha=0.1)

    fill_between([data_min, 0.0, -data_min],
        [-data_min, -data_min, -data_min],
        [-data_min, 0.0, -data_min],
        color="grey",
        lw=0.0,
        alpha=0.1)


    xlabel("observed slip rate (mm a⁻¹)")
    ylabel("modeled slip rate (mm a⁻¹)")
    title("extension")

    tight_layout()

    return fig
end


function plot_resid_dists(obs, pred; obs_err=nothing, pred_err=nothing, norm=false,
    histogram=false, cum_hist=false, pct_lines = (0.1, 0.25, 0.5, 0.75, 0.9), abs=false,
    trim_axes=false, trim_val=0.99)

    resids = obs .- pred

    if abs
        resids = abs.(resids)
    end

    sort_idx = sortperm(resids)

    x_arr = collect(1:length(obs))
    
    if (isnothing(obs_err) & isnothing(pred_err))
        errs = nothing
    elseif !isnothing(obs_err) & isnothing(pred_err)
        errs = obs_err
    elseif !isnothing(pred_err) & isnothing(obs_err)
        errs = pred_err
    else
        errs = sqrt.(obs_err.^2 .+ pred_err.^2)
    end
    
    if norm == false
        resid_plot = resids[sort_idx]
        if !isnothing(errs)
            errs_plot = errs[sort_idx]
        else
            errs_plot = errs
        end
    else
        resid_plot = (resids ./ abs.(obs))[sort_idx]
        if !isnothing(errs)
            errs_plot = (errs./ abs.(obs))[sort_idx]
        else
            errs_plot = errs
        end
    end

    mean_resids = sum(resid_plot) / length(resids)
    sd_resids = std(resid_plot)
    println("mean resids: $mean_resids")
    println("std resids: $sd_resids")

    if (histogram == false) & (cum_hist == false)
        #errorbar(x_arr, resid_plot, yerr=errs_plot, fmt=",", elinewidth=0.3)
        errorbar(resid_plot, x_arr, xerr=errs_plot, fmt=",", elinewidth=0.3)
    elseif (histogram == false) & (cum_hist == true)
        n_bins = Int(ceil(length(resid_plot) / 5))
        hist(resid_plot, bins=n_bins, cumulative=true)
    elseif (histogram == true) & (cum_hist == false)
        n_bins = Int(ceil(length(resid_plot) / 5))
        hist(resid_plot, bins=n_bins)
    elseif (histogram == true) & (cum_hist == true)
        n_bins = Int(ceil(length(resid_plot) / 5))
        hist(resid_plot, bins=n_bins, density=false)
        ax2 = gca().twinx()
        plot(resid_plot, collect(1:length(obs))./length(obs),
             color="C1", lw=0.5)
        for pct in pct_lines
            axhline(pct, color="grey", lw=0.25)
        end
    end

end


function plot_geol_resids(dex_geol_obs, dex_geol_err, dex_geol_pred, ext_geol_obs, 
    ext_geol_err, ext_geol_pred)

    fig = figure(figsize=(7,7))
    
    subplot(3,2,1)
    plot_resid_dists(dex_geol_obs, dex_geol_pred; obs_err=dex_geol_err)
    xlabel("residuals (mm/yr)")
    title("strike-slip residuals")

    subplot(3,2,3)
    plot_resid_dists(dex_geol_obs, dex_geol_pred; obs_err=dex_geol_err, norm=true)
    xlabel("relative residuals")
    title("strike-slip residuals\n(normalized to observation rate)")
    
    subplot(3,2,5)
    plot_resid_dists(dex_geol_obs, dex_geol_pred; obs_err=dex_geol_err, histogram=true)
    xlabel("residuals (mm/yr)")
    title("strike-slip residuals")

    subplot(3,2,2)
    plot_resid_dists(ext_geol_obs, ext_geol_pred; obs_err=ext_geol_err)
    xlabel("residuals (mm/yr)")
    title("extension residuals")

    subplot(3,2,4)
    plot_resid_dists(ext_geol_obs, ext_geol_pred; obs_err=ext_geol_err, norm=true)
    xlabel("relative residuals")
    title("extension residuals\n(normalized to observation rate)")

    return fig
end


function plot_gnss_component_resids(e_obs, e_err, e_pred, n_obs, 
    n_err, n_pred)

    fig = figure(figsize=(7,7))
    
    subplot(2,2,1)
    plot_resid_dists(n_obs, n_pred; obs_err=n_err)
    ylabel("residuals (mm/yr)")
    title("Vn residuals")

    subplot(2,2,3)
    plot_resid_dists(n_obs, n_pred; obs_err=n_err, norm=true)
    ylabel("relative residuals")
    title("Vn residuals\n(normalized to observation rate)")

    subplot(2,2,2)
    plot_resid_dists(e_obs, e_pred; obs_err=e_err)
    ylabel("residuals (mm/yr)")
    title("Ve residuals")

    subplot(2,2,4)
    plot_resid_dists(e_obs, e_pred; obs_err=e_err, norm=true)
    ylabel("relative residuals")
    title("Ve residuals\n(normalized to observation rate)")

    return fig
end


function plot_gnss_vector_resids(e_obs, e_err, e_pred, n_obs, 
    n_err, n_pred)

    e_misfit = e_obs .- e_pred
    n_misfit = n_obs .- n_pred

    pred_mags = sqrt.(e_obs.^2 .+ n_obs.^2)

    misfits, misfit_errs = Oiler.Geom.rotate_xy_vec_to_magnitude(
        e_misfit, n_misfit; x_err=e_err, y_err=n_err
    )

    println(misfits[1:5])
    println(misfit_errs[1:5])

    fig = figure()#figsize=(7,7))

    subplot(2,1,1)
    plot_resid_dists(misfits, zeros(size(misfits)); obs_err=misfit_errs,
        histogram=true, cum_hist=true)
        xlabel("GNSS misfit magnitude (mm/yr)")

    subplot(2,1,2)
    plot_resid_dists((misfits ./ pred_mags), zeros(size(misfits)); obs_err=misfit_errs,
        histogram=true, cum_hist=true, norm=false)
        xlabel("normalized GNSS misfit magnitude (mm/yr)")

    return fig
end


function plot_residuals(results)

    vels = results["data"]["gnss_results"]
    vel_fig = plot_gnss_component_resids(vels.obs_ve, vels.obs_ee, vels.pred_ve,
    vels.obs_vn, vels.obs_en, vels.pred_vn)

    # obs_ve, pred_ve
    #vel_df.obs_ve, vel_df.obs_ee, vel_df.obs_ee)

    #plot_geol_resids()
    #return vel_fig

    vel_mag_fig = plot_gnss_vector_resids(vels.obs_ve, vels.obs_ee, vels.pred_ve,
    vels.obs_vn, vels.obs_en, vels.pred_vn)

    return vel_fig, vel_mag_fig
    #return vel_mag_fig
end


end # module

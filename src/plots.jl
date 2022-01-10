module Plots

using PyPlot
using DataFrames, DataFramesMeta

using Oiler



function plot_results_map(results, vel_groups, faults, tris=[])
    vel_df = Oiler.ResultsAnalysis.get_gnss_results(results, vel_groups)

    fig = figure(figsize=(14, 14))

    for fault in faults
        plot(fault.trace[:,1], fault.trace[:,2], "k-", lw=0.5)
    end

    if length(tris) > 0
        cm = get_cmap(:viridis)
    
        tri_rates = [get_tri_total_rate(tri, results) for tri in tris]
        tri_rate_min = minimum(tri_rates)
        tri_rate_max = maximum(tri_rates)
    
        #clim(tri_rate_min, tri_rate_max)
    
        for tri in tris
            plot_tri(tri, results; vmin = tri_rate_min, vmax = tri_rate_max, cm = cm)
        end
        #colorbar()
    end

    quiver(vel_df.lon, vel_df.lat,
           vel_df.obs_ve, vel_df.obs_vn,
           color="b", scale=300, label="obs GNSS")
    quiver(vel_df.lon, vel_df.lat,
           vel_df.pred_ve, vel_df.pred_vn,
           color="r", scale=300, label="pred GNSS")

    legend()
    axis("equal")
    
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
        # fault_idx = parse(Int, rate.name)
        fault_idx = rate.name
        fault_row = @where(fault_df, :fid .== fault_idx)[1,:]
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
    
    
    inc = map(!, iszero.(geol_slip_rate_df[!,:include]))
    
    function check_missing(val)
        if val == ""
            return false
        elseif val == 0.
            return false
        elseif ismissing(val) | isnothing(val)
            return false
        else
            return true
        end
    end
    
    dex_inc = [check_missing(d) for d in geol_slip_rate_df[!,:dextral_rate]]
    ext_inc = [check_missing(d) for d in geol_slip_rate_df[!,:extension_rate]]
    
    dex_geol_obs = geol_slip_rate_df[!,:dextral_rate][dex_inc .* inc]
    # dex_geol_err = parse.(Float64, geol_slip_rate_df[!,:dextral_err][dex_inc .* inc])
    dex_geol_err = geol_slip_rate_df[!,:dextral_err][dex_inc .* inc];
    dex_geol_pred = [p[1] for p in pred_geol_slip_rates[dex_inc[inc]]];
    dex_geol_pred_err = [p[3] for p in pred_geol_slip_rates[dex_inc[inc]]];
    
    ext_geol_obs = geol_slip_rate_df[!,:extension_rate][ext_inc .* inc];
    ext_geol_err = geol_slip_rate_df[!,:extension_err][ext_inc .* inc];
    ext_geol_pred =  [p[2] for p in pred_geol_slip_rates[ext_inc[inc]]];
    ext_geol_pred_err = [p[4] for p in pred_geol_slip_rates[ext_inc[inc]]];
    
    fig = figure(figsize=(5, 9));
    # suptitle("Observed vs. Modeled Quaternary Slip Rates")
    
    subplot(2, 1, 1);
    data_max = maximum([maximum(dex_geol_obs), maximum(dex_geol_pred)]);
    data_min = minimum([minimum(dex_geol_obs), minimum(dex_geol_pred)]);
    
    plot([data_min, data_max], [data_min, data_max], "C1--", lw=0.5);
    
    axis("equal");
    errorbar(dex_geol_obs, dex_geol_pred, xerr=dex_geol_err, yerr=dex_geol_pred_err,
             fmt=",", elinewidth=0.3);
    
    autoscale(false);
    
    fill_between([data_min, 0., -data_min], 
                 [data_min, 0., data_min], 
                 [data_min, data_min, data_min],
                 color="grey",
                 lw=0.,
                 alpha=0.1);
    
    fill_between([data_min, 0., -data_min], 
                 [-data_min, -data_min, -data_min],
                 [-data_min, 0., -data_min], 
                 color="grey",
                 lw=0.,
                 alpha=0.1);
    
    xlabel("observed");
    ylabel("modeled");
    title("dextral");
    subplot(2, 1, 2);
    
    data_max = maximum([maximum(ext_geol_obs), maximum(ext_geol_pred)]);
    data_min = minimum([minimum(ext_geol_obs), minimum(ext_geol_pred)]);
    plot([data_min, data_max], [data_min, data_max], "C1--", lw=0.5);
    
    axis("equal")
    errorbar(ext_geol_obs, ext_geol_pred, xerr=ext_geol_err, yerr=ext_geol_pred_err, 
             fmt=",", elinewidth=0.3);
    
    autoscale(false);
    
    fill_between([data_min, 0., -data_min], 
                 [data_min, 0., data_min], 
                 [data_min, data_min, data_min],
                 color="grey",
                 lw=0.,
                 alpha=0.1);
    
    fill_between([data_min, 0., -data_min], 
                 [-data_min, -data_min, -data_min],
                 [-data_min, 0., -data_min], 
                 color="grey",
                 lw=0.,
                 alpha=0.1);
    
    
    xlabel("observed");
    ylabel("modeled");
    title("extension");
    
    tight_layout();
    
    return fig
end


end # module

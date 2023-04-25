module Elastic

using Logging
using ThreadsX

import Base.Threads.@spawn
import Base.Threads.@threads

using ..Oiler
using ..Oiler: Fault, VelocityVectorSphere, fault_oblique_merc, get_gnss_vels,
    get_coords_from_vel_array, az_to_angle, okada



export fault_to_okada

"""
    fault_to_okada(fault::Fault, sx1::Float64, sy1::Float64, 
    sx2::Float64, sy2::Float64)

Formats the attributes of a `Fault` into a form suitable for calculating
dislocations using Okada's equations.

# Arguments
- fault: fault that will be used to calculate dislocations
- sx1

"""
function fault_to_okada(fault::Fault, sx1::Float64, sy1::Float64,
    sx2::Float64, sy2::Float64)

    usd = fault.usd * 1000.0
    lsd = fault.lsd * 1000.0

    strike = atan(sy1 - sy2, sx1 - sx2) + pi
    # strike = az_to_angle(fault.strike) 
    delta = deg2rad(fault.dip)
    L = sqrt((sx2 - sx1)^2 + (sy2 - sy1)^2)
    W = (lsd - usd) / sin(delta)

    # calculate fault segment anchor and other buried point
    ofx = sx1 + lsd / tan(delta) * sin(strike)
    ofy = sy1 - lsd / tan(delta) * cos(strike)

    Dict("strike" => strike, "delta" => delta, "L" => L, "W" => W, "lsd" => lsd,
        "ofx" => ofx, "ofy" => ofy)
end


function calc_okada_locking_effects_per_fault(fault::Fault, lons, lats;
    elastic_floor=1e-4, check_nans=false, fill_nans=true, nan_fill_val=0.0,
    ss=1.0, ds=1.0, ts=1.0)

    # project fault and velocity coordinates
    n_gnss = length(lons)
    xp, yp = Oiler.Faults.fault_proj(fault, lons, lats)
    xg, yg = xp[1:n_gnss], yp[1:n_gnss] # gnss
    sx1, sy1, sx2, sy2 = xp[end-1], yp[end-1], xp[end], yp[end] # fault 

    # format okada
    d = fault_to_okada(fault, sx1, sy1, sx2, sy2)

    # calc okada partials
    elastic_partials = distribute_partials(okada(d, ss, ts, ds,
        xg, yg; floor=elastic_floor))

    if check_nans == true
        for (i, ep) in enumerate(elastic_partials)
            if any(isnan, ep)
                @warn "NaN"
                @warn i
                @warn fault
                bad_lon, bad_lat = lons[i], lats[i]
                @warn "$bad_lon, $bad_lat"
                if fill_nans == true
                    map!(x -> isnan(x) ? nan_fill_val : x, ep, ep)
                    if any(isnan, ep)
                        @warn "Couldn't replace NaN"
                    end
                end
            end
        end
    end

    # build rotation and transformation matrices
    Pf = Oiler.Faults.build_velocity_projection_matrix(fault.strike, fault.dip)

    PvGbf = Oiler.BlockRotations.build_PvGb_vel(
        Oiler.fault_to_vel(fault)
    )

    PfPvGbf = Pf * PvGbf

    [part * PfPvGbf for part in elastic_partials]
end


function calc_locking_effects_segmented_fault(fault::Fault, lons, lats;
    elastic_floor=1e-4, check_nans=false)
    # may have some problems w/ dip dir for highly curved faults
    trace = fault.trace
    simp_trace = Oiler.Geom.simplify_polyline(trace, 0.01)
    # simp_trace = trace

    parts = []
    for i = 1:size(simp_trace, 1)-1
        seg_trace = simp_trace[i:i+1, :]
        part = calc_okada_locking_effects_per_fault(Oiler.Fault(trace=seg_trace,
                dip=fault.dip,
                dip_dir=fault.dip_dir,
                lsd=fault.lsd,
                usd=fault.usd,
                check_trace=false,
                #usd = maximum([fault.usd; 1.0])
            ),
            lons, lats; elastic_floor=elastic_floor,
            check_nans=check_nans)
        push!(parts, part)
    end
    sum(parts)
end


function calc_any_fault_locking_effects(fault::Fault, lons, lats;
    elastic_floor=1e-4, check_nans=false, dip_threshold=0.0)

    # tri locking calcs are slower and may be a bit more inaccurate;
    # for now they are turned off.

    if fault.dip >= dip_threshold
        return calc_locking_effects_segmented_fault(fault, lons, lats;
            elastic_floor=elastic_floor, check_nans=check_nans)
    else
        return calc_tri_locking_effects_per_fault(fault, lons, lats;
            elastic_floor=elastic_floor, check_nans=check_nans)
    end
end

function get_tri_slip_vec(ds, ss, fault_strike, tri_strike)
    angle_diff = Oiler.Geom.angle_difference(fault_strike, tri_strike, return_abs=false)
    rad_diff = deg2rad(angle_diff)
    new_ds, new_ss = Oiler.Geom.rotate_velocity(ds, ss, rad_diff)
end

function calc_tri_locking_effects_per_fault(fault::Fault, lons, lats;
    ss=1.0, ds=1.0, ts=0.0, elastic_floor=1e-4, check_nans=false, fill_nans=true,
    nan_fill_val=0.0)

    trace = fault.trace
    simp_trace = Oiler.Geom.simplify_polyline(trace, 0.2)
    simp_fault = Fault(trace=simp_trace, dip=fault.dip, lsd=fault.lsd,
        dip_dir=fault.dip_dir, check_trace=false)

    fault_tris = get_fault_tris(simp_fault)

    ds_slip_vecs = zeros(length(fault_tris))
    ss_slip_vecs = zeros(length(fault_tris))

    for (i, tri) in enumerate(fault_tris)
        tri_strike = Oiler.Tris.get_tri_strike_dip(tri)[1]
        ds_slip_vecs[i], ss_slip_vecs[i] = get_tri_slip_vec(ds, ss, simp_fault.strike, tri_strike)
    end

    due = zeros(length(lons))
    dun = zeros(length(lons))
    duv = zeros(length(lons))
    sue = zeros(length(lons))
    sun = zeros(length(lons))
    suv = zeros(length(lons))
    txx = zeros(length(lons))

    for (i, tri) in enumerate(fault_tris)
        de, dn, dv, se, sn, sv = calc_tri_slip(tri, lons, lats;
            ss=ss_slip_vecs[i], ds=ds_slip_vecs[i], floor=elastic_floor)
        due .+= de
        dun .+= dn
        duv .+= dv
        sue .+= se
        sun .+= sn
        suv .+= sv
    end

    elastic_partials = distribute_partials((sue, sun, suv,
        due, dun, duv, txx, txx, txx))

    if check_nans == true
        for (i, ep) in enumerate(elastic_partials)
            if any(isnan, ep)
                @warn "NaN"
                @warn i
                @warn fault
                bad_lon, bad_lat = lons[i], lats[i]
                @warn "$bad_lon, $bad_lat"
                if fill_nans == true
                    map!(x -> isnan(x) ? nan_fill_val : x, ep, ep)
                    if any(isnan, ep)
                        @warn "Couldn't replace NaN"
                    end
                end
            end
        end
    end

    # build rotation and transformation matrices
    Pf = Oiler.Faults.build_velocity_projection_matrix(fault.strike, fault.dip)

    PvGbf = Oiler.BlockRotations.build_PvGb_vel(
        Oiler.fault_to_vel(fault)
    )

    PfPvGbf = Pf * PvGbf

    [part * PfPvGbf for part in elastic_partials]
end


function get_fault_tris(fault::Fault)
    lower_trace = Oiler.Faults.project_fault_trace(fault)

    tris = []
    for i = 1:size(fault.trace, 1)-1
        t1 = [fault.trace[i, 1]; fault.trace[i, 2]; 0.0]
        t2 = [fault.trace[i+1, 1]; fault.trace[i+1, 2]; 0.0]
        b1 = [lower_trace[i, 1]; lower_trace[i, 2]; -fault.lsd]
        b2 = [lower_trace[i+1, 1]; lower_trace[i+1, 2]; -fault.lsd]
        tri1 = Oiler.Tris.Tri(p1=t1, p2=t2, p3=b1)
        tri2 = Oiler.Tris.Tri(p1=t2, p2=b1, p3=b2)
        push!(tris, tri1)
        push!(tris, tri2)
    end

    tris
end



function distribute_partials(partials::NTuple{9,Array{Float64}})
    # (es, ens, )
    n_gnss = length(partials[1])
    [make_partials_matrix(partials, i) for i = 1:n_gnss]
end


function make_partials_matrix(partials, i::Integer)
    # es ed et
    # ns nd nt
    # us ud ut 

    # the vertical velocities are not returned, as this will cause problems
    # when vertical velocities from GNSS are not used, as is typical.

    [partials[1][i] partials[4][i] partials[7][i]
        partials[2][i] partials[5][i] partials[8][i]
        # partials[3][i] partials[6][i] partials[9][i]]
        0.0 0.0 0.0]
end


function calc_locking_effects(faults, vel_groups; elastic_floor=1e-4,
    check_nans=false)

    gnss_vels = get_gnss_vels(vel_groups)
    gnss_lons = [vel["vel"].lon for vel in gnss_vels]
    gnss_lats = [vel["vel"].lat for vel in gnss_vels]
    gnss_idxs = [vel["idx"] for vel in gnss_vels]

    vg_keys = sort(collect(Tuple(keys(vel_groups))))
    fault_groups = Oiler.Utils.group_faults(faults, vg_keys)
    locking_partial_groups = Dict()

    # calculate locking effects from each fault at each site
    # locking effects for faults in each vel_group sum
    #for vg in vg_keys
    @threads for vg in vg_keys
        if haskey(fault_groups, vg)
            locking_partial_groups[vg] = sum([
                #calc_locking_effects_segmented_fault(fault, gnss_lons, gnss_lats,
                calc_any_fault_locking_effects(fault, gnss_lons, gnss_lats,
                    elastic_floor=elastic_floor, check_nans=check_nans)
                for fault in fault_groups[vg]
            ])
        else
        end
    end

    # make a new dictionary with keys as the indices of the sub-array in
    # the model matrix
    locking_partials = Dict()
    for (i, vg) in enumerate(vg_keys)
        if haskey(locking_partial_groups, vg)
            col_id = 3 * (i - 1) + 1
            col_idx = col_id:col_id+2
            for (j, row_idx) in enumerate(gnss_idxs)
                locking_partials[(row_idx[1], col_idx)] = locking_partial_groups[vg][j]
            end
        end
    end
    # @info "done calculating locking partials"
    locking_partials
end

"""
    calc_tri_effects

Calculates unit strike and dip slip on all GNSS velocities for all tris.
Returns an 3Mx2T matrix, where M is the number of sites and T is the number
of tris. Only horizontal components are returned.
"""
function calc_tri_effects(tris, gnss_lons, gnss_lats; elastic_floor=1e-4)

    tri_gnss_partials = hcat(ThreadsX.collect(
    #tri_gnss_partials = hcat(collect(
        arrange_tri_partials(
            calc_tri_effects_single_tri(tri, gnss_lons, gnss_lats;
                elastic_floor=elastic_floor)...)
        for tri in tris
    )...)
end


function calc_tri_effects_single_tri(tri, lons, lats; elastic_floor=1e-4)

    ss = 1.# e-3 # mm
    ds = 1.# e-3 # mm
    ts = 0.0 # no tensile slip on tris

    calc_tri_slip(tri, lons, lats; ss=ss, ds=ds,
        ts=ts, floor=elastic_floor)
end


function calc_tri_slip(tri, lons, lats; ss=0.0, ds=0.0, ts=0.0, floor=1e-4)

    # project coordinates of tri and sites
    x_proj, y_proj = Oiler.Tris.tri_proj(tri, lons, lats)
    x_gnss, y_gnss = x_proj[1:length(lons)], y_proj[1:length(lats)]
    z_gnss = zeros(length(x_gnss))

    # make projected tri
    tri_x_1, tri_y_1 = x_proj[end-2], y_proj[end-2]
    tri_x_2, tri_y_2 = x_proj[end-1], y_proj[end-1]
    tri_x_3, tri_y_3 = x_proj[end], y_proj[end]

    tri_z_1 = tri.p1[3] * 1000.0
    tri_z_2 = tri.p2[3] * 1000.0
    tri_z_3 = tri.p3[3] * 1000.0

    tri_p1 = [tri_x_1 tri_y_1 tri_z_1]
    tri_p2 = [tri_x_2 tri_y_2 tri_z_2]
    tri_p3 = [tri_x_3 tri_y_3 tri_z_3]

    # dip slip component
    due, dun, duv = Oiler.TD.TDdispHS(x_gnss, y_gnss, z_gnss, tri_p1, tri_p2,
        tri_p3, 0.0, ds)
    due[abs.(due).<abs(ds * floor)] .= 0.0
    dun[abs.(dun).<abs(ds * floor)] .= 0.0
    duv[abs.(duv).<abs(ds * floor)] .= 0.0

    # strike_component
    sue, sun, suv = Oiler.TD.TDdispHS(x_gnss, y_gnss, z_gnss, tri_p1, tri_p2,
        tri_p3, ss, 0.0)

    sue[abs.(sue).<abs(ss * floor)] .= 0.0
    sun[abs.(sun).<abs(ss * floor)] .= 0.0
    suv[abs.(suv).<abs(ss * floor)] .= 0.0

    due', dun', duv', sue', sun', suv'
end


function arrange_tri_partials(due, dun, duv, sue, sun, suv; uv_zero=true)

    # one 3x2 matrix for each site, stacked vertically
    # due sue
    # dun sun
    # duv suv

    n_vels = length(due)

    if uv_zero == true
        duv = zeros(size(duv))
        suv = zeros(size(suv))
    end

    out = zeros(n_vels * 3, 2)

    for i = 1:n_vels
        ind_v = i * 3
        ind_n = ind_v - 1
        ind_e = ind_v - 2

        out[ind_e, 1] = due[i]
        out[ind_n, 1] = dun[i]
        out[ind_v, 1] = duv[i]
        out[ind_e, 2] = sue[i]
        out[ind_n, 2] = sun[i]
        out[ind_v, 2] = suv[i]
    end
    out
end

end # module

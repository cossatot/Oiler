module Elastic

using Logging
using ThreadsX

import Base.Threads.@spawn
import Base.Threads.@threads

using ..Oiler
using ..Oiler: Fault, okada



const PartialMatrix = Matrix{Float64}
const PartialStack = Vector{PartialMatrix}
const LockingPartialKey = Tuple{UnitRange{Int},UnitRange{Int}}


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


function build_fault_transform_matrix(fault::Fault)
    Pf = Oiler.Faults.build_velocity_projection_matrix(fault.strike, fault.dip)
    vlon, vlat = Oiler.Faults.get_midpoint(fault.trace)
    PvGbf = Oiler.BlockRotations.build_PvGb_deg(vlon, vlat)

    Pf * PvGbf
end


function transform_partial_stack(elastic_partials::PartialStack, fault::Fault)
    PfPvGbf = build_fault_transform_matrix(fault)
    transformed = Vector{Matrix{Float64}}(undef, length(elastic_partials))

    @inbounds for i in eachindex(elastic_partials)
        transformed[i] = elastic_partials[i] * PfPvGbf
    end

    transformed
end


function accumulate_partial_stacks!(dest::PartialStack, src::PartialStack)
    @inbounds for i in eachindex(dest, src)
        dest[i] .+= src[i]
    end

    dest
end


function collect_gnss_locking_inputs(vel_groups, vg_keys)
    gnss_lons = Float64[]
    gnss_lats = Float64[]
    gnss_rows = UnitRange{Int}[]
    row_set_num = 0

    for key in vg_keys
        for vel in vel_groups[key]
            row_set_num += 1
            if vel.vel_type == "GNSS"
                row_idx = 3 * (row_set_num - 1) + 1
                push!(gnss_lons, vel.lon)
                push!(gnss_lats, vel.lat)
                push!(gnss_rows, row_idx:row_idx+2)
            end
        end
    end

    gnss_lons, gnss_lats, gnss_rows
end


function group_faults_by_key(faults, vel_group_keys)
    fault_groups = Dict{Tuple{String,String},Vector{Fault}}()

    for fault in faults
        for key in vel_group_keys
            if fault.fw in key && fault.hw in key
                push!(get!(() -> Fault[], fault_groups, key), fault)
            end
        end
    end

    fault_groups
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
                    #map!(x -> isnan(x) ? nan_fill_val : x, ep, ep)
                    #Oiler.Utils.replacenan!(ep, replacement=nan_fill_val)
                    replace!(ep, NaN=>nan_fill_val)
                    if any(isnan, ep)
                        @warn "Couldn't replace NaN"
                    end
                end
            end
        end
    end

    transform_partial_stack(elastic_partials, fault)
end


function calc_locking_effects_segmented_fault(fault::Fault, lons, lats;
    elastic_floor=1e-4, check_nans=false)
    # may have some problems w/ dip dir for highly curved faults
    trace = fault.trace
    simp_trace = Oiler.Geom.simplify_polyline(trace, 0.01)
    # simp_trace = trace

    n_segments = size(simp_trace, 1) - 1
    if n_segments <= 0
        return [zeros(3, 3) for _ in eachindex(lons)]
    end

    parts = calc_okada_locking_effects_per_fault(Oiler.Fault(trace=simp_trace[1:2, :],
            dip=fault.dip,
            dip_dir=fault.dip_dir,
            lsd=fault.lsd,
            usd=fault.usd,
            check_trace=false,
        ),
        lons, lats; elastic_floor=elastic_floor,
        check_nans=check_nans)

    for i = 2:n_segments
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
        accumulate_partial_stacks!(parts, part)
    end

    parts
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
                    #map!(x -> isnan(x) ? nan_fill_val : x, ep, ep)
                    replace!(ep, NaN=>nan_fill_val)
                    if any(isnan, ep)
                        @warn "Couldn't replace NaN"
                    end
                end
            end
        end
    end

    transform_partial_stack(elastic_partials, fault)
end


function get_fault_tris(fault::Fault)
    lower_trace = Oiler.Faults.project_fault_trace(fault)

    n_segs = size(fault.trace, 1) - 1
    tris = Vector{Oiler.Tris.Tri}(undef, 2 * n_segs)

    for i = 1:n_segs
        t1 = [fault.trace[i, 1]; fault.trace[i, 2]; 0.0]
        t2 = [fault.trace[i+1, 1]; fault.trace[i+1, 2]; 0.0]
        b1 = [lower_trace[i, 1]; lower_trace[i, 2]; -fault.lsd]
        b2 = [lower_trace[i+1, 1]; lower_trace[i+1, 2]; -fault.lsd]
        tri1 = Oiler.Tris.Tri(p1=t1, p2=t2, p3=b1)
        tri2 = Oiler.Tris.Tri(p1=t2, p2=b1, p3=b2)
        tri_idx = 2 * i - 1
        tris[tri_idx] = tri1
        tris[tri_idx+1] = tri2
    end

    tris
end



function distribute_partials(partials::NTuple{9,Vector{Float64}})
    # (es, ens, )
    n_gnss = length(partials[1])
    out = Vector{Matrix{Float64}}(undef, n_gnss)

    @inbounds for i = 1:n_gnss
        out[i] = make_partials_matrix(partials, i)
    end

    out
end


function make_partials_matrix(partials, i::Integer)
    # es ed et
    # ns nd nt
    # us ud ut 

    # the vertical velocities are not returned, as this will cause problems
    # when vertical velocities from GNSS are not used, as is typical.

    out = Matrix{Float64}(undef, 3, 3)
    out[1, 1] = partials[1][i]
    out[1, 2] = partials[4][i]
    out[1, 3] = partials[7][i]
    out[2, 1] = partials[2][i]
    out[2, 2] = partials[5][i]
    out[2, 3] = partials[8][i]
    out[3, 1] = 0.0
    out[3, 2] = 0.0
    out[3, 3] = 0.0

    out
end

function orient_fault_to_key(fault, key)
    if (fault.hw, fault.fw) == key
        return fault
    elseif (fault.fw, fault.hw) == key
        rev_trace = reverse(fault.trace, dims=1)
        new_dip_dir = Oiler.Faults.get_closest_dir(
            Oiler.Faults.get_trace_dip_trend_rhr(rev_trace)
        )
        return Oiler.Fault(
            trace=rev_trace,
            dip=fault.dip,
            dip_dir=new_dip_dir,
            extension_rate=fault.extension_rate,
            extension_err=fault.extension_err,
            dextral_rate=fault.dextral_rate,
            dextral_err=fault.dextral_err,
            cde=fault.cde,
            lsd=fault.lsd,
            usd=fault.usd,
            name=fault.name,
            hw=fault.fw,
            fw=fault.hw,
            fid=fault.fid,
            check_trace=false,
        )
    else
        @warn "Fault $(fault.fid) does not match velocity group $key"
        return fault
    end
end

function calc_locking_effects(faults, vel_groups; elastic_floor=1e-4,
    check_nans=false)

    @info "\tprepping vels"
    vg_keys = sort(collect(Tuple(keys(vel_groups))))
    gnss_lons, gnss_lats, gnss_rows = collect_gnss_locking_inputs(vel_groups, vg_keys)

    @info "\tsorting keys and grouping faults"
    fault_groups = group_faults_by_key(faults, vg_keys)
    locking_partial_groups = Dict{Tuple{String,String},PartialStack}()


    # calculate locking effects from each fault at each site
    # locking effects for faults in each vel_group sum
    @info "\tcalculating locking per fault group (parallel?)"
    #for vg in vg_keys
    n_threads = Threads.nthreads()
    @info "$n_threads threads"
    #@threads for vg in vg_keys
    for vg in vg_keys
        if haskey(fault_groups, vg)
            oriented_faults = [orient_fault_to_key(fault, vg) for fault in fault_groups[vg]]
            n_faults = length(oriented_faults)
            if n_faults > 0
                group_partials = calc_any_fault_locking_effects(oriented_faults[1],
                    gnss_lons, gnss_lats, elastic_floor=elastic_floor,
                    check_nans=check_nans)
                for fault_idx = 2:n_faults
                    fault_partials = calc_any_fault_locking_effects(oriented_faults[fault_idx],
                        gnss_lons, gnss_lats,
                        elastic_floor=elastic_floor, check_nans=check_nans)
                    accumulate_partial_stacks!(group_partials, fault_partials)
                end
                locking_partial_groups[vg] = group_partials
            end
        else
        end
    end

    # make a new dictionary with keys as the indices of the sub-array in
    # the model matrix
    @info "\t\tmaking partials dicts"
    locking_partials = Dict{LockingPartialKey,Matrix{Float64}}()
    for (i, vg) in enumerate(vg_keys)
        if haskey(locking_partial_groups, vg)
            col_id = 3 * (i - 1) + 1
            col_idx = col_id:col_id+2
            for (j, row_idx) in enumerate(gnss_rows)
                locking_partials[(row_idx, col_idx)] = locking_partial_groups[vg][j]
            end
        end
    end
    @info "\tdone calculating locking partials"
    locking_partials
end

"""
    calc_tri_effects

Calculates unit strike and dip slip on all GNSS velocities for all tris.
Returns an 3Mx2T matrix, where M is the number of sites and T is the number
of tris. Only horizontal components are returned.
"""
function calc_tri_effects(tris, gnss_lons, gnss_lats; elastic_floor=1e-4)

    if isempty(tris)
        return zeros(length(gnss_lons) * 3, 0)
    end

    tri_gnss_partials = hcat(ThreadsX.collect(
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

    due, dun, duv, sue, sun, suv
end


function arrange_tri_partials(due, dun, duv, sue, sun, suv; uv_zero=true)

    # one 3x2 matrix for each site, stacked vertically
    # due sue
    # dun sun
    # duv suv

    n_vels = length(due)

    out = zeros(n_vels * 3, 2)

    for i = 1:n_vels
        ind_v = i * 3
        ind_n = ind_v - 1
        ind_e = ind_v - 2

        out[ind_e, 1] = due[i]
        out[ind_n, 1] = dun[i]
        out[ind_v, 1] = uv_zero ? 0.0 : duv[i]
        out[ind_e, 2] = sue[i]
        out[ind_n, 2] = sun[i]
        out[ind_v, 2] = uv_zero ? 0.0 : suv[i]
    end
    out
end


end # module

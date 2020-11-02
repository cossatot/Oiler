module Elastic

import Base.Threads.@spawn
import Base.Threads.@threads

using ..Oiler
using ..Oiler: Fault, VelocityVectorSphere, fault_oblique_merc, get_gnss_vels,
    get_coords_from_vel_array, az_to_angle, okada


using Logging

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

    usd = fault.usd * 1000.
    lsd = fault.lsd * 1000.

    strike = atan(sy1 - sy2, sx1 - sx2) + pi 
    # strike = az_to_angle(fault.strike) 
    delta = deg2rad(fault.dip)
    L = sqrt((sx2 - sx1)^2 + (sy2 - sy1)^2)
    W = (lsd - usd) / sin(delta)

    # calculate fault segment anchor and other buried point
    ofx =  sx1 + lsd / tan(delta) * sin(strike)
    ofy =  sy1 - lsd / tan(delta) * cos(strike)

    Dict("strike" => strike, "delta" => delta, "L" => L, "W" => W, "lsd" => lsd,
         "ofx" => ofx, "ofy" => ofy)
end


function calc_locking_effects_per_fault(fault::Fault, lons, lats)

    # project fault and velocity coordinates
    n_gnss = length(lons)
    xp, yp = fault_oblique_merc(fault, lons, lats)
    xg, yg = xp[1:n_gnss], yp[1:n_gnss] # gnss
    sx1, sy1, sx2, sy2 = xp[end - 1], yp[end - 1], xp[end], yp[end] # fault 

    # TODO: add distance-based filtering HERE to maintain some sparsity in
    # the results to preserve RAM

    # format Okada
    D = fault_to_okada(fault, sx1, sy1, sx2, sy2)

    # calc Okada partials
    elastic_partials = distribute_partials(okada(D, 1., 1., 1.,   
        xg, yg))

    # build rotation and transformation matrices
    Pf = Oiler.Faults.build_velocity_projection_matrix(fault.strike, fault.dip)

    PvGbf = Oiler.BlockRotations.build_PvGb_vel(
        Oiler.fault_to_vel(fault)
    )

    PfPvGbf = Pf * PvGbf

    [part * PfPvGbf for part in elastic_partials]
end


function calc_locking_effects_segmented_fault(fault::Fault, lons, lats)
    # may have some problems w/ dip dir for highly curved faults
    trace = fault.trace
    simp_trace = Oiler.Geom.simplify_polyline(trace, 0.2)
    # simp_trace = trace

    parts = []
    for i in 1:size(simp_trace, 1) - 1
        seg_trace = simp_trace[i:i + 1,:]
        part = calc_locking_effects_per_fault(Oiler.Fault(trace=seg_trace, 
                dip=fault.dip, 
                dip_dir=fault.dip_dir,
                lsd=fault.lsd, usd=fault.usd), 
            lons, lats)
        push!(parts, part)
    end
    sum(parts)
end


function distribute_partials(partials::NTuple{9,Array{Float64}})
    n_gnss = length(partials[1])
    [make_partials_matrix(partials, i) for i in 1:n_gnss]
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
     0. 0. 0.]
end


function calc_locking_effects(faults, vel_groups)

    gnss_vels = get_gnss_vels(vel_groups)
    gnss_lons = [vel["vel"].lon for vel in gnss_vels]
    gnss_lats = [vel["vel"].lat for vel in gnss_vels]
    gnss_idxs = [vel["idx"] for vel in gnss_vels]

    vg_keys = sort(collect(Tuple(keys(vel_groups))))
    fault_groups = Oiler.Utils.group_faults(faults, vg_keys)
    locking_partial_groups = Dict()

    # calculate locking effects from each fault at each site
    # locking effects for faults in each vel_group sum
    @threads for vg in vg_keys
        if haskey(fault_groups, vg)
            locking_partial_groups[vg] = sum([
                calc_locking_effects_segmented_fault(fault, gnss_lons, gnss_lats)
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
            col_idx = col_id:col_id + 2
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
Returns an Mx2T matrix, where M is the number of sites and T is the number
of tris. Only horizontal components are returned.
"""
function calc_tri_effects(tris, gnss_lons, gnss_lats)

    tri_gnss_partials = hcat( collect(
        [arrange_tri_partials(
            calc_tri_effects_single_tri(tri, gnss_lons, gnss_lats)...)
         for tri in tris]
    )...)

end


function calc_tri_effects_single_tri(tri, lons, lats)

    # project coordinates of tri and sites
    x_proj, y_proj = Oiler.Tris.tri_merc(tri, lons, lats)
    x_gnss, y_gnss = x_proj[1:length(lons)], y_proj[1:length(lats)]
    z_gnss = zeros(length(x_gnss))

    # make projected tri
    tri_x_1, tri_y_1 = x_proj[end - 2], y_proj[end - 2]
    tri_x_2, tri_y_2 = x_proj[end - 1], y_proj[end - 1]
    tri_x_3, tri_y_3 = x_proj[end], y_proj[end]

    tri_z_1 = tri.p1[3] * 1000.
    tri_z_2 = tri.p2[3] * 1000.
    tri_z_3 = tri.p3[3] * 1000.

    tri_p1 = [tri_x_1 tri_y_1 tri_z_1]
    tri_p2 = [tri_x_2 tri_y_2 tri_z_2]
    tri_p3 = [tri_x_3 tri_y_3 tri_z_3]

    ss_slip = 1.# e-3 # mm
    ds_slip = 1.# e-3 # mm
    ts_slip = 0. # no tensile slip on tris

    # dip slip component
    due, dun, duv = Oiler.TD.TDdispHS(x_gnss, y_gnss, z_gnss, tri_p1, tri_p2, 
                                     tri_p3, 0., ds_slip)
    
    # strike_slip_component
    sue, sun, suv = Oiler.TD.TDdispHS(x_gnss, y_gnss, z_gnss, tri_p1, tri_p2, 
                                     tri_p3, ss_slip, 0.)
    due', dun', duv', sue', sun', suv'
end


function arrange_tri_partials(due, dun, duv, sue, sun, suv; uv_zero=true)
    n_vels = length(due)

    if uv_zero == true
        duv = zeros(size(duv))
        suv = zeros(size(duv))
    end

    out = zeros(n_vels * 3, 2)

    for i in 1:n_vels
        ind_v = i * 3
        ind_n = ind_v - 1
        ind_e = ind_v - 2

        out[ind_e,1] = due[i]
        out[ind_n,1] = dun[i]
        out[ind_v,1] = duv[i]
        out[ind_e,2] = sue[i]
        out[ind_n,2] = sun[i]
        out[ind_v,2] = suv[i]
    end
    out
end
end # module

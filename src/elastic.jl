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
    ofxe = sx2 + lsd / tan(delta) * sin(strike)
    ofye = sy2 - lsd / tan(delta) * cos(strike)

    # calculate fault segment anchor and other buried point (top relative)
    tfx =  sx1 + usd / tan(delta) * sin(strike)
    tfy =  sy1 - usd / tan(delta) * cos(strike)
    tfxe = sx2 + usd / tan(delta) * sin(strike)
    tfye = sy2 - usd / tan(delta) * cos(strike)

    Dict("strike" => strike, "delta" => delta, "L" => L, "W" => W, "lsd" => lsd,
         "ofx" => ofx, "ofy" => ofy, "ofxe" => ofxe, "ofye" => ofye,
         "tfx" => tfx, "tfy" => tfy, "tfxe" => tfxe, "tfye" => tfye)
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
    #simp_trace = trace

    parts = []
    for i in 1:size(simp_trace,1) - 1
        seg_trace = simp_trace[i:i+1,:]
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
    [partials[1][i] partials[4][i] partials[7][i]
     partials[2][i] partials[5][i] partials[8][i]
     #partials[3][i] partials[6][i] partials[9][i]]
     0. 0. 0.]
end


function calc_locking_effects(faults, vel_groups)

    #@info "calculating locking partials"

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

end # module

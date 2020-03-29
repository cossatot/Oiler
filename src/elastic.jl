module Elastic

using ..Oiler: Fault, VelocityVectorSphere, fault_oblique_merc, get_gnss_vels,
    get_coords_from_vel_array, az_to_angle

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

    # strike = az_to_angle(fault.strike) 
    strike = atan(sy1 - sy2, sx1 - sx2) + pi 
    delta = deg2rad(fault.dip)
    L = sqrt((sx2 - sx1)^2 + (sy2 - sy1)^2)
    W = (fault.lsd - fault.usd) / sin(delta)

    # calculate fault segment anchor and other buried point
    ofx = sx1 + fault.lsd / tan(delta) * sin(strike)
    ofy = sy1 - fault.lsd / tan(delta) * cos(strike)
    ofxe = sx2 + fault.lsd / tan(delta) * sin(strike)
    ofye = sy2 - fault.lsd / tan(delta) * cos(strike)

    # calculate fault segment anchor and other buried point (top relative)
    tfx = sx1 + fault.usd / tan(delta) * sin(strike)
    tfy = sy1 - fault.usd / tan(delta) * cos(strike)
    tfxe = sx2 + fault.usd / tan(delta) * sin(strike)
    tfye = sy2 - fault.usd / tan(delta) * cos(strike)

    # just to keep everything lined up nicely below
    lsd = fault.lsd

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

    # format Okada
    D = fault_to_okada(fault, sx1, sy1, sx2, sy2)

    # calc Okada
    

end

function calc_locking_effects(faults, vel_groups)
    gnss_vels = get_gnss_vels(vel_groups)
    gnss_lons, gnss_lats = get_coords_from_vel_array(gnss_vels)

    fault_lock_displs = sum([calc_locking_effects_per_fault(fault, lons, lats)
        for fault in faults])

    # for each row/station
    #     calc Pf, project
    #     sparsify
    
end

function add_effects_to_PvGb()
end

end # module
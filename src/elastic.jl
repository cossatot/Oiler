module Elastic

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


"""
    okada_wrapper(fault::Dict, strike_disp::Float64, dip_disp::Float64, 
    tensile_disp::Float64, xs::Array{Float64}, ys::Array{Float64}, 
    Pr::Float64 = 0.25)

Translates the outputs from the Okada's dislocation results into a more
appropriate form for Oiler. The ENU format is translated into the NED reference
frame that Oiler uses, the results are put into a matrix, and a weird
backwards thing is fixed.
"""
function okada_wrapper(fault::Dict, strike_disp::Float64, dip_disp::Float64, 
    tensile_disp::Float64, xs::Array{Float64}, ys::Array{Float64}, 
    Pr::Float64 = 0.25)

    (ves, vns, vus, ved, vnd, vud, vet, vnt, vut) = okada(fault, strike_disp,
        dip_disp, tensile_disp, xs, ys, Pr
    )


end

function calc_locking_effects_per_fault(fault::Fault, lons, lats)

    # project fault and velocity coordinates
    n_gnss = length(lons)
    xp, yp = fault_oblique_merc(fault, lons, lats)
    xg, yg = xp[1:n_gnss], yp[1:n_gnss] # gnss
    sx1, sy1, sx2, sy2 = xp[end - 1], yp[end - 1], xp[end], yp[end] # fault 

    # format Okada
    D = fault_to_okada(fault, sx1, sy1, sx2, sy2)

    # calc Okada partials
    elastic_partials = distribute_partials(okada(D, -1., -1., -1., xg, yg))

    # build rotation and transformation matrices
    Pf = Oiler.Faults.build_velocity_projection_matrix(fault.strike, fault.dip)

    PvGbf = Oiler.BlockRotations.build_PvGb_vel(
        Oiler.fault_to_vel(fault)
    )

    PfPvGbf = Pf * PvGbf

    [part * PfPvGbf for part in elastic_partials]
end


function distribute_partials(partials::NTuple{9,Array{Float64}})
    n_gnss = length(partials[1])
    [make_partials_matrix(partials, i) for i in 1:n_gnss]
end


function make_partials_matrix(partials, i::Integer)
    # [partials[1][i] partials[2][i] partials[3][i]
    # partials[4][i] partials[5][i] partials[6][i]
    # partials[7][i] partials[8][i] partials[9][i]]
    [partials[1][i] partials[4][i] partials[7][i]
     partials[2][i] partials[5][i] partials[8][i]
     partials[3][i] partials[6][i] partials[9][i]]
end


function calc_locking_effects(faults, vel_groups, sparse_tol::Float64 = 1e-5)
    gnss_vels = get_gnss_vels(vel_groups)
    gnss_lons, gnss_lats = get_coords_from_vel_array(gnss_vels)

    # fault_lock_displs = sum([calc_locking_effects_per_fault(fault, gnss_lons, gnss_lats)
    #    for fault in faults])

    # for each row/station
    #     calc Pf, project
    #     sparsify
    
end

function add_effects_to_PvGb()
end

end # module
module Okada

export okada

using ..Oiler:rotate_xy_vec

"""
    okada(fault::Dict, strike_disp::Float64, dip_disp::Float64, 
        tensile_disp::Float64, xs::Array{Float64}, ys::Array{Float64}, 
        Pr::Float64 = 0.25)

Calculates the slip in an elastic halfspace due to a rectangular
dislocation.

Ported from 'okada_partials.m' by Meade and Loveless.

# Arguments
- `fault`: A dictionary of fault parameters that has been output by
    `fault_to_okada`.
- `strike_disp`: The strike-slip displacement.
- `dip_dips`: The dip-slip displacment.
- `tens_disp`: The tensile displacement.
- `pr`: Poisson ratio

# Returns
- `uxstot`: slip in `x` direction due to strike slip
- `uystot`: slip in `y` direction due to strike slip
- `uzstot`: slip in `z` direction due to strike slip
- `uxdtot`: slip in `x` direction due to dip slip
- `uydtot`: slip in `y` direction due to dip slip
- `uzdtot`: slip in `z` direction due to dip slip
- `uxttot`: slip in `x` direction due to tensile slip
- `uyttot`: slip in `y` direction due to tensile slip
- `uzttot`: slip in `z` direction due to tensile slip
"""
function okada(fault::Dict, strike_disp::Float64, dip_disp::Float64, 
    tensile_disp::Float64, xs::Array{Float64}, ys::Array{Float64}; 
    Pr::Float64=0.25, floor=0.)
    
    # unpack dict and convert variable names, all inherited from Blocks
    xf = fault["ofx"]
    yf = fault["ofy"]
    strike = fault["strike"]
    d = fault["lsd"]
    delta = fault["delta"]
    L = fault["L"]
    W = fault["W"]

    U1 = strike_disp
    U2 = dip_disp
    U3 = tensile_disp

    # constants
    tol = 1e-4
    alpha = -2. * Pr + 1.

    # station locations relative to the fault anchor
    xt = xs .- xf
    yt = ys .- yf

    # rotate stations to remove strike
    alpha_rot = -strike
    xr, yr = rotate_xy_vec(xt, yt, alpha_rot)

    # calculate some values that are frequently needed
    sindd = sin(delta)
    cosdd = cos(delta)
    twopi = 2. * pi

    nx = length(xr)

    # uxtot = zeros(nx, 1)
    # uytot = zeros(nx, 1)
    # uztot = zeros(nx, 1)

    # find displacements at each station
    x = xr
    y = yr
    p = y .* cosdd .+ d .* sindd
    q = repeat(y .* sindd .- d .* cosdd, 1, 4)
    zi = [x x x .- L x .- L]
    eta = [p p .- W p p .- W]
    ybar = eta .* cosdd + q .* sindd
    dbar = eta .* sindd - q .* cosdd
    R = sqrt.(zi.^2 .+ eta.^2 .+ q.^2)
    X = sqrt.(zi.^2 + q.^2)

    # Calculate some more commonly used values.
    # These are introduced to reduce repetitive
    # calculations.
    Reta = R + eta
    Rzi = R + zi
    Rdbar = R + dbar
    qdivR = q ./ R
    phi = atan.(zi .* eta ./ q ./ R)

    if abs(cosdd) >= tol
        I5 = (alpha .* 2 ./ cosdd .* atan.((eta .* (X .+ q .* cosdd) 
                .+ (X .* (R + X) * sindd)) ./ (zi .* (R .+ X) .* cosdd)) )

        I4 = alpha ./ cosdd * (log.(Rdbar) .- sindd .* log.(Reta))

        I3 = (alpha .* (1. ./ cosdd .* ybar ./ Rdbar .- log.(Reta) ) .+ sindd 
              ./ cosdd .* I4)

        I2 = alpha .* (-log.(Reta)) .- I3
        I1 = alpha .* (-1. ./ cosdd .* zi ./ Rdbar) .- sindd ./ cosdd .* I5

    else
        I5 = -alpha .* (zi .* sindd) ./ Rdbar
        I4 = -alpha .* q ./ Rdbar
        I3 = alpha ./ 2 .* (eta ./ Rdbar .+ ybar .* q ./ Rdbar.^2 .- log.(Reta))
        I2 = alpha .* (-log.(Reta)) .- I3
        I1 = -alpha ./ 2. .* (zi .* q) ./ Rdbar.^2 
    end # if

    uxs = -U1 ./ twopi .* (zi .* qdivR ./ Reta .+ phi + I1 .* sindd)
    uxd = -U2 ./ twopi .* (qdivR .- I3 .* sindd .* cosdd)
    uxt =  U3 ./ twopi .* (q .* qdivR ./ Reta .- I3 .* sindd.^2)

    uys = -U1 ./ twopi .* (ybar .* qdivR ./ Reta .+ q .* cosdd ./ Reta + I2 .* sindd)
    uyd = -U2 ./ twopi .* (ybar .* qdivR ./ Rzi  .+ cosdd .* phi .- I1 .* sindd .* cosdd)
    uyt =  (U3 ./ twopi .* (-dbar .* qdivR ./ Rzi .- sindd .* (zi .* qdivR ./ Reta
            .- phi) .- I1 .* sindd.^2))

    
    uzs = -U1 ./ twopi .* (dbar .* qdivR ./ Reta + q .* sindd ./ Reta + I4 .* sindd)
    uzd = -U2 ./ twopi .* (dbar .* qdivR ./ Rzi + sindd .* phi .- I5 .* sindd .* cosdd)
    uzt =  (U3 ./ twopi .* (ybar .* qdivR ./ Rzi .+ cosdd .* (zi .* qdivR ./ Reta .-
            phi) .- I5 .* sindd.^2))
         
    uxstot = uxs[:, 1] - uxs[:, 2] - uxs[:, 3] + uxs[:, 4]
    uxdtot = uxd[:, 1] - uxd[:, 2] - uxd[:, 3] + uxd[:, 4]
    uxttot = uxt[:, 1] - uxt[:, 2] - uxt[:, 3] + uxt[:, 4]
    uystot = uys[:, 1] - uys[:, 2] - uys[:, 3] + uys[:, 4]
    uydtot = uyd[:, 1] - uyd[:, 2] - uyd[:, 3] + uyd[:, 4]
    uyttot = uyt[:, 1] - uyt[:, 2] - uyt[:, 3] + uyt[:, 4]
    uzstot = uzs[:, 1] - uzs[:, 2] - uzs[:, 3] + uzs[:, 4]
    uzdtot = uzd[:, 1] - uzd[:, 2] - uzd[:, 3] + uzd[:, 4]
    uzttot = uzt[:, 1] - uzt[:, 2] - uzt[:, 3] + uzt[:, 4]

    # rotate the station displacements back to include the effects of the strike
    uxstot, uystot  = rotate_xy_vec(uxstot, uystot, -alpha_rot)
    uxdtot, uydtot  = rotate_xy_vec(uxdtot, uydtot, -alpha_rot)
    uxttot, uyttot  = rotate_xy_vec(uxttot, uyttot, -alpha_rot)

    # replace very small values with 0. to preserve sparsity
    tot_disp = sqrt(U1^2 + U2^2 + U3^2)
    for vv in (uxstot, uystot, uzstot, uxdtot, 
               uydtot, uzdtot, uxttot, uyttot, uzttot)
        vv[abs.(vv) .< abs(floor * tot_disp)] .= 0.
    end

    uxstot, uystot, uzstot, uxdtot, uydtot, uzdtot, uxttot, uyttot, uzttot
end



end # module
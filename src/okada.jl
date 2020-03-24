module Okada


"""
    okada

Calculates the slip in an elastic halfspace due to a rectangular
dislocation.

Ported from 'okada_partials.m' by Meade and Loveless.
"""
function okada(xf, yf, strike, d, delta, L, W, U1, U2, U3, xs, ys, Pr)

    # constants
    tol = 1e-4
    alpha = -2. * Pr + 1.

    # station locations relative to the fault anchor
    xt = xs - xf
    yt = ys - yf

    # rotate stations to remove strike
    alpha_rot = -strike
    xr, yr = rotate_xy_vec(xt, yt, alpha_rot)

    # calculate some valuees that are fequently needed
    sindd = sin(delta)
    cosdd = cos(delta)
    twopi = 2. * pi

    nx = length(xr)

    #uxtot = zeros(nx, 1)
    #uytot = zeros(nx, 1)
    #uztot = zeros(nx, 1)

    # find displacements at each station
    x = xr
    y = yr
    p = y .* cosdd + d .* sindd
    q = repeat(y .* sindd - d .* cosdd, 1, 4)
    zi = [x x x-L x-L]
    eta = [p p-W p p-W]
    ybar = eta .* cosdd + q .* sindd
    dbar = eta .* sindd - q .* cosdd
    R = sqrt(zi.^2 + eta.^2 + q.^2)
    X = sqrt(zi.^2 + q.^2)

    # Calculate some more commonly used values.
    # These are introduced to reduce repetitive
    # calculations.
    Reta = R + eta
    Rzi = R + zi
    Rdbar = R + dbar
    qdivR = q ./ R
    phi = atan(zi .* eta ./ q ./ R)

    if abs(cosdd) >= tol
        I5 = (alpha * 2 ./ cosdd * atan((eta .* (X + q .* cosdd) 
                + (X .* (R + X) * sindd)) ./ (zi .* (R+X) .* cosdd) ) )

        I4 = alpha ./ cosdd * (log(Rdbar) - sindd .* log(Reta))

        I3 = (alpha .* (1. ./ cosdd .* ybar ./ Rdbar - log(Reta) ) + sindd 
              ./ cosdd .* I4)

        I2 = alpha .* (-log(Reta)) - I3
        I1 = alpha .* (-1. ./ cosdd .* zi ./ Rdbar) - sindd ./ cosdd .* I5

    else
        I5 = -alpha .* (zi .* sindd) ./ Rdbar
        I4 = -alpha .* q ./ Rdbar
        I3 = alpha ./ 2 .* (eta ./ Rdbar + ybar .* q ./ Rdbar.^2 - log(Reta))
        I2 = alpha .* (-log(Reta) - I3
        I1 = -alpha/2. .* (zi .* q) ./ Rdbar.^2 
    end #if

    uxs = -U1 ./ twopi .* (zi .* qdivR ./ Reta + phi + I1 .* sindd)
    uxd = -U2 ./ twopi .* (qdivR - I3 .* sindd .* cosdd)
    uxt =  U3 ./ twopi .* (q .* qdivR ./ Reta - I3 .* sindd.^2)

    uys = -U1 ./ twopi .* (ybar .* qdivR ./ Reta + q .* cosdd ./ Reta + I2 .* sindd)
    uyd = -U2 ./ twopi .* (ybar .* qdivR ./ Rzi)  + cosd .* phi - I1 .* sind .* cosd)
    uyt =  (U3 ./ twopi .* (-dbar .* qdivR ./ Rzi - sind .* (zi .* qdivR ./ Reta
            - phi) - I1 .* sind.^2))

    
    uzs = -U1 ./ twopi .* (dbar .* qdivR ./ Reta + q .* sind ./ Reta + I4 .* sindd);
    uzd = -U2 ./ twopi .* (dbar .* qdivR ./ Rzi + sind .* phi - I5 .* sind .* cosdd)
    uzt =  U3 ./ twopi .* (ybar .* qdivR ./ Rzi + cosd .* zi .* qdivR ./ Reta -
    phi) - I5 .* sind.^2)
         
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

end


function rotate_xy_vec(x, y, alpha_rot)
    xp = cosd(alpha_rot) .* x - sind(alpha_rot) .* y
    yp = sind(alpha_rot) .* x + cosd(alpha_rot) .* y
    return xp, yp
end

end #module
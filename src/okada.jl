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

    uxtot = zeros(nx, 1)
    uytot = zeros(nx, 1)
    uztot = zeros(nx, 1)

    # find displacements at each station
    x = xr
    y = yr


    


end


function rotate_xy_vec(xt, yt, alpha_rot)
end

end
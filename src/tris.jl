module Tris

export Tri

using Setfield
using LinearAlgebra

import Proj: Transformation#Projection, transform

using ..Oiler


struct Tri
    p1::Array{Float64,1}
    p2::Array{Float64,1}
    p3::Array{Float64,1}
    dip_slip_rate::Float64
    dip_slip_err::Float64
    strike_slip_rate::Float64
    strike_slip_err::Float64
    cds::Float64
    dip_locking_frac::Float64
    strike_locking_frac::Float64
    name::String
    hw::String
    fw::String
end


function Tri(; 
    p1::Array{Float64,1},
    p2::Array{Float64,1},
    p3::Array{Float64,1},
    dip_slip_rate::Float64=0.,
    dip_slip_err::Float64=0.,
    strike_slip_rate::Float64=0.,
    strike_slip_err::Float64=0.,
    cds::Float64=0.,
    dip_locking_frac::Float64=0.,
    strike_locking_frac::Float64=0.,
    name::String="",
    hw::String="",
    fw::String=""
    )

    for z in [p1[3] p2[3] p3[3]]
        if z > 0.0
            warn_msg = "Tri $name has Z coordinate $z > 0 (above HS surface)"
            @warn warn_msg
        end
    end

    if (p1 == p2) || (p1 == p3) || (p2 == p3)
        warn_msg = "Tri $name has duplicated coordinates"
        @warn warn_msg
    end

    if Oiler.Geom.check_winding_order([p1, p2, p3, p1]) == 1
        # @warn "reversing tri"
        p1, p2, p3 = (p3, p2, p1)
    end

    Tri(p1, p2, p3, dip_slip_rate, dip_slip_err, strike_slip_rate, 
        strike_slip_err, cds, dip_locking_frac, strike_locking_frac, name, hw, 
        fw)
end


"""
    get_tri_center()
Returns the centroid of the tri.
"""
function get_tri_center(tri::Oiler.Tris.Tri)
    #lon1 = tri.p1[1]
    #lat1 = tri.p1[2]
    #lon2 = tri.p2[1]
    #lat2 = tri.p2[2]
    #lon3 = tri.p3[1]
    #lat3 = tri.p3[2]
    #z1 = tri.p1[3]
    #z2 = tri.p2[3]
    #z3 = tri.p3[3]

    #lons = [lon1 lon2 lon3]
    #lats = [lat1 lat2 lat3]
    #zs = [z1 z2 z3]
    
    #ps1, ps2 = get_tri_strike_line(tri.p1, tri.p2, tri.p3)
    
    #wgs84 = "+proj=longlat +datum=WGS84 +nodefs"
    #omerc = Oiler.Geom.get_oblique_merc(ps1[1], ps1[2], ps2[1], ps2[2])
    #trans = Transformation(wgs84, omerc; always_xy=true)

    #xy = [trans(lon, lats[i]) for (i, lon) in enumerate(lons)]
    #xs = [c[1] for c in xy]
    #ys = [c[2] for c in xy]

    #xs, ys = tri_merc(tri, lons, lats)

    p1c = Oiler.Geom.point_sphere_to_cart(tri.p1)
    p2c = Oiler.Geom.point_sphere_to_cart(tri.p2)
    p3c = Oiler.Geom.point_sphere_to_cart(tri.p3)

    xs = [p1c[1] p2c[1] p3c[1]]
    ys = [p1c[2] p2c[2] p3c[2]]
    zs = [p1c[3] p2c[3] p3c[3]]


    xc = sum(xs) / 3.
    yc = sum(ys) / 3.
    zc = sum(zs) / 3.

    lon, lat, depth = Oiler.Geom.point_cart_to_sphere(xc, yc, zc)

    #untrans = Transformation(omerc, wgs84; always_xy=true)
    #ll = untrans(xc, yc)

    #[ll[1], ll[2], zc]
end


"""
    tri_merc()

Performs a Mercator projection localized on the surface projection of the
first two points of the tri.

"""
function tri_merc(tri, lons, lats)

    lon1 = tri.p1[1]
    lat1 = tri.p1[2]
    lon2 = tri.p2[1]
    lat2 = tri.p2[2]
    lon3 = tri.p3[1]
    lat3 = tri.p3[2]
    z1 = tri.p1[3]
    z2 = tri.p2[3]
    z3 = tri.p3[3]

    lons_w_tri_pts = lons[:]
    lats_w_tri_pts = lats[:]
    append!(lons_w_tri_pts, [lon1 lon2 lon3])
    append!(lats_w_tri_pts, [lat1 lat2 lat3])

    ps1, ps2 = get_tri_strike_line(tri.p1, tri.p2, tri.p3)

    xp, yp = Oiler.Geom.oblique_merc(lons_w_tri_pts, lats_w_tri_pts, ps1[1],
        ps1[2], ps2[1], ps2[2])
end


function tri_azimuthal(tri, lons, lats)
    lon1 = tri.p1[1]
    lat1 = tri.p1[2]
    lon2 = tri.p2[1]
    lat2 = tri.p2[2]
    lon3 = tri.p3[1]
    lat3 = tri.p3[2]
    z1 = tri.p1[3]
    z2 = tri.p2[3]
    z3 = tri.p3[3]

    lons_w_tri_pts = lons[:]
    lats_w_tri_pts = lats[:]
    append!(lons_w_tri_pts, [lon1 lon2 lon3])
    append!(lats_w_tri_pts, [lat1 lat2 lat3])

    ps1, ps2 = get_tri_strike_line(tri.p1, tri.p2, tri.p3)

    xp, yp = Oiler.Geom.azimuthal_equidistant_proj(lons_w_tri_pts, lats_w_tri_pts, ps1[1],
        ps1[2], ps2[1], ps2[2])
end


function tri_proj(tri, lons, lats; proj=:azimuthal)
    if proj == :mercator
        xp, yp = tri_merc(tri, lons, lats)
    elseif proj == :azimuthal
        xp, yp = tri_azimuthal(tri, lons, lats)
    end

    xp, yp
end


"""
    get_tri_strike_line()

Returns two [lat, lon, depth] points from three triangular vertices
such that the two points are on different sides of the triangle but
are at the same depth; a line between them would be a line of strike.

If two of the triangle's vertices are at equal depths, these points
are returned. Otherwise, the middle depth vertex and a point on the
side opposite are returned.
"""
function get_tri_strike_line(p1, p2, p3)
    if p1[3] == p2[3]
        strike_line = (p1, p2)
    elseif p1[3] == p3[3]
        strike_line = (p1, p3)
    elseif p2[3] == p3[3]
        strike_line = (p2, p3)
    else
        # sort by z coordinate (depth is negative)
        p_low, p_med, p_high = sort([p1, p2, p3], by=x -> x[end])

        # find point on other leg at point of equal depth to p_med
        # this point is at a horizontal fraction of the great circle distance
        # between the points that is proportional to the depth difference
        # between the bottom and mid points
        Z = p_low[3] - p_high[3]
        z = p_low[3] - p_med[3]
        L = Oiler.Geom.gc_distance(p_high[1], p_high[2], p_low[1], p_low[2])
        l = L * z / Z

        p_mid_other_leg = Oiler.Geom.sample_polyline(
            vcat(vec(p_low)', vec(p_high)'), [l])[1]
    
        strike_line = (p_med, [p_mid_other_leg[1], p_mid_other_leg[2], p_med[3]])
    end

    strike_line
end



function get_tri_strike_line(tri::Tri)
    get_tri_strike_line(tri.p1, tri.p2, tri.p3)
end


function get_tri_strike_dip(tri::Tri)
    proj_x, proj_y = tri_merc(tri, [0.], [0.])
    xp1 = [proj_x[1] proj_y[1] tri.p1[3] * 1000]
    xp2 = [proj_x[2] proj_y[2] tri.p2[3] * 1000]
    xp3 = [proj_x[3] proj_y[3] tri.p3[3] * 1000]

    strike, dip = Oiler.Geom.strike_dip_from_3_pts(xp1, xp2, xp3)
end


function get_tri_rate_from_pole(tri::Tri, pole)
    center = Oiler.Tris.get_tri_center(tri)
    pred_vel = Oiler.BlockRotations.predict_block_vel(center[1], center[2], pole)

    strike, dip = get_tri_strike_dip(tri)

    v_rl, v_ex = Oiler.Faults.ve_vn_to_fault_slip_rate(pred_vel.ve, pred_vel.vn, 
                                                       strike)
    
    if (pred_vel.ee == 0.) & (pred_vel.en == 0.)
        e_rl, e_ex, cde = 0., 0., 0.
    else
        e_rl, e_ex, cde = Oiler.Faults.ee_en_to_fault_slip_rate_err(pred_vel.ee, 
                                                  pred_vel.en, strike; 
                                                  cen=pred_vel.cen)
    end
    strike_slip_rate = v_rl
    strike_slip_err = e_rl
    dip_slip_rate = -v_ex
    dip_slip_err = e_ex
    cds = cde # maybe backwards?

    dip_slip_rate, dip_slip_err, strike_slip_rate, strike_slip_err, cds
end

    
"""
    
"""
function tri_centroid_distance(tri1::Oiler.Tris.Tri, tri2::Oiler.Tris.Tri)
    c1 = get_tri_center(tri1)
    c2 = get_tri_center(tri2)
    
    c1_cart = Oiler.Geom.point_sphere_to_cart(c1)
    c2_cart = Oiler.Geom.point_sphere_to_cart(c2)
    
    sqrt(sum((c1_cart - c2_cart).^2))
end


function point_colocation(p1, p2; horiz_decimal=4, vert_decimal=1)
    (round(p1[1], digits=horiz_decimal) == round(p2[1], digits=horiz_decimal)) & 
    (round(p1[2], digits=horiz_decimal) == round(p2[2], digits=horiz_decimal)) & 
    (round(p1[3], digits=vert_decimal) == round(p2[3],  digits=vert_decimal))
end


"""
    check_tri_adjacence()

Checks to see if two tris are adjacent, i.e. sharing a single vertex or an
edge.

# Arguments
    - `tri1`: An Oiler.Tris.Tri instance
    - `tri2`: An Oiler.Tris.Tri instance
    - `n_common_pts`: Number of points in common to define adjacence. Should be 
      `1` (a singl shared vertex means tris are adjacent) or `2` (a shared
       edge is required for tri adjacence). Defaults to `2`.
    - `self_adjacence`: A boolean to determine whether a tri can be adjacent
       to itself, or another tri with the same geometry. Defaults to `false`.

# Returns
    - bool, `true` or `false`
"""
function check_tri_adjacence(tri1, tri2; n_common_pts::Int=2, 
                             self_adjacence=false,
                             horiz_decimal=5, vert_decimal=1)
    if (n_common_pts != 1) & (n_common_pts != 2)
        throw(ArgumentError("$n_common_pts must be 1 or 2"))
    end

    common_pts = 0
    for pp in [tri1.p1, tri1.p2, tri1.p3]
        for cc in [tri2.p1, tri2.p2, tri2.p3]
            if point_colocation(pp, cc)
                common_pts += 1
            end
        end
    end

    if !self_adjacence
        if common_pts == 3
            return false
        end
    end

    if common_pts < n_common_pts
        return false
    else
        return true
    end
end


function get_tri_adjacence_dict(tris; n_common_pts=2, self_adjacence=false)
    tri_adj_dict = Dict()
    for t1 in tris
        tri_adj_dict[t1.name] = []
        for t2 in tris
            if Oiler.Tris.check_tri_adjacence(t1, t2; 
                                              n_common_pts=n_common_pts, 
                                              self_adjacence=self_adjacence)
                push!(tri_adj_dict[t1.name], t2.name)
            end
        end
    end
    tri_adj_dict
end

# useless without results
# function get_tri_total_rate(tri::Oiler.Tris.Tri)
#    ds = results["tri_slip_rates"][tri.name]["dip_slip"]
#    ss = results["tri_slip_rates"][tri.name]["strike_slip"]
#    total_rate = sqrt(ds^2 + ss^2)
# end


function get_tri_hanging_wall(tri::Tri, block_df; epsg=4326, set=true)
    hw = Oiler.IO.get_block_idx_for_point(get_tri_center(tri), block_df, 
                                          epsg=epsg)
    if set
        tri = @set tri.hw = hw
        return tri
    else
        return hw
    end
end




end # module

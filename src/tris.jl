module Tris

export Tri

using ..Oiler


struct Tri
    p1::Array{Float64,1}
    p2::Array{Float64,1}
    p3::Array{Float64,1}
    dip_slip_rate::Float64
    strike_slip_rate::Float64
    name::String
end


function Tri(; 
    p1::Array{Float64,1},
    p2::Array{Float64,1},
    p3::Array{Float64,1},
    dip_slip_rate::Float64=0.,
    strike_slip_rate::Float64=0.,
    name::String=""
    )

    if Oiler.Geom.check_winding_order([p1, p2, p3]) == 1
        println("reversing tri")
        p1, p2, p3 = (p3, p2, p1)
    end

    Tri(p1, p2, p3, dip_slip_rate, strike_slip_rate, name)
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

    if lat1 == 0.
        lat1 += 1e-4
    end

    ps1, ps2 = get_tri_strike_line(tri.p1, tri.p2, tri.p3)

    xp, yp = Oiler.Geom.oblique_merc(lons_w_tri_pts, lats_w_tri_pts, ps1[1], 
                                ps1[2], ps2[1], ps2[2])
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
        # println(p1)
        # println(p2)
        # println(p3)
        # sort by z coordinate (depth is positive)
        p_high, p_med, p_low = sort([p1, p2, p3], by=x -> x[end])

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
        println(p_mid_other_leg)
    
        strike_line = (p_med, [p_mid_other_leg[1], p_mid_other_leg[2], p_med[3]])
    end
    strike_line
end


end # module
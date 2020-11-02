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

    # use tri corners for projection center
    xp, yp = Oiler.Geom.oblique_merc(lons_w_tri_pts, lats_w_tri_pts, lon1, lat1, 
                                      lon2, lat2)
end


end # module
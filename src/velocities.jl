module Velocities

export VelocityVectorSphere

using Parameters



"""
    VelocityVectorSph(20., )
Velocity vector in spherical (lon, lat) coordinates, with velocities in
mm/yr.

ve, vn, vu are east, north and up velocities, and ee, en, and eu are the 
1-sigma uncertainties.
"""
@with_kw struct VelocityVectorSphere
    lond::Float64
    latd::Float64
    ve::Float64
    vn::Float64
    vd::Float64 = 0.
    ee::Float64 = 0.
    en::Float64 = 0.
    ed::Float64 = 0.
    fix::String = ""
    mov::String = ""
    name::String = ""
end

function reverse(vel::VelocityVectorSphere)
    VelocityVectorSphere(lond = vel.lond,
        latd = vel.latd,
        ve = -vel.ve,
        vn = -vel.vn,
        vd = -vel.vd,
        ee = vel.ee,
        en = vel.en,
        ed = vel.ed,
        fix = vel.mov,
        mov = vel.fix,
        name = name)
end


function Base.:-(vel::VelocityVectorSphere)
    reverse(vel)
end

end
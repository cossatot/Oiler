module RotationPoles

export PoleCart, PoleSphere, add_poles, add_poles, subtract_poles,
    pole_cart_to_sphere, pole_sphere_to_cart

using LinearAlgebra

using Parameters


"""
Euler Pole (rotation vector) in Cartesian coordinates

# Arguments
- `x::Float64`: *x* coordinate in the Cartesian system
- `y::Float64`: *y* coordinate in the Cartesian system
- `z::Float64`: *z* coordinate in the Cartesian system
- `fix::String`: Name of fixed block. Defaults to "".
- `mov::String`: Name of moving (movative) block.  Defaults to "".
"""
@with_kw struct PoleCart
    x::Float64
    y::Float64
    z::Float64
    ex::Float64 = 0.
    ey::Float64 = 0.
    ez::Float64 = 0.
    fix::String = ""
    mov::String = ""
end


"""
Euler Pole (rotation vector) in spherical coordinates.

# Arguments
- `lon::Float64`: Longitude coordinate in the spherical system (in degrees)
- `lat::Float64`: Latitude coordinate in the spherical system (in degrees)
- `rotrate::Float64`: Rotation rate (in degree / time unit, unspecified)
- `fix::String`: Name of fixed block. Defaults to "".
- `mov::String`: Name of moving (movative) block.  Defaults to "".
"""
@with_kw struct PoleSphere
    lon::Float64
    lat::Float64
    rotrate::Float64
    elon::Float64 = 0.
    elat::Float64 = 0.
    erotrate::Float64 = 0.
    fix::String = ""
    mov::String = ""
end


function Base.:-(pole::PoleCart)
    PoleCart(pole; x=-pole.x, y=-pole.y, z=-pole.z)
end

function Base.:-(pole::PoleSphere)
    pc = pole_sphere_to_cart(pole)
    pole_cart_to_sphere(-pc)
end


function add_poles(poles::PoleCart...)
    xn = sum(pole.x for pole in poles)
    yn = sum(pole.y for pole in poles)
    zn = sum(pole.z for pole in poles)

    PoleCart(x=xn, y=yn, z=zn, fix=poles[1].fix, 
                  mov=poles[end].mov)
end


function add_poles(poles::Array{PoleCart})
    xn = sum(pole.x for pole in poles)
    yn = sum(pole.y for pole in poles)
    zn = sum(pole.z for pole in poles)

    PoleCart(x=xn, y=yn, z=zn, fix=poles[1].fix, 
                  mov=poles[end].mov)
end


function add_poles(poles::PoleSphere...)
    cart_poles = [pole_sphere_to_cart(pole) for pole in poles]

    cart_pole_sum = add_poles(cart_poles...)
    pole_cart_to_sphere(cart_pole_sum)
end


function Base.:+(pole1::PoleCart, pole2::PoleCart)
    add_poles(pole1, pole2)
end

function Base.:+(pole1::PoleSphere, pole2::PoleSphere)
    add_poles(pole1, pole2)
end

function subtract_poles(pole1::PoleCart, pole2::PoleCart)
    xx = pole1.x - pole2.x
    yy = pole1.y - pole2.y
    zz = pole1.z - pole2.z

    PoleCart(x=xx, y=yy, z=zz, fix=pole1.fix, mov=pole2.fix)
end

function subtract_poles(pole1::PoleSphere, pole2::PoleSphere)
    pole1c = pole_sphere_to_cart(pole1)
    pole2c = pole_sphere_to_cart(pole2)

    pole_diff = subtract_poles(pole1c, pole2c)

    pole_cart_to_sphere(pole_diff)
end


function Base.:-(pole1::PoleCart, pole2::PoleCart)
    subtract_poles(pole1, pole2)
end

function Base.:-(pole1::PoleSphere, pole2::PoleSphere)
    subtract_poles(pole1, pole2)
end


function err_cart_to_sphere_bm(x::Float64, y::Float64, z::Float64, 
    ex::Float64, ey::Float64, ez::Float64)
    # Derivation modified slightly from 'OmegaSigToEulerSig.m' in the Blocks
    # package (Meade and Loveless, 2009 BSSA)

    # short-circuit if errors are not present to prevent numerical inaccuracies
    # or NANs when the values are zero
    if all(x -> x == 0., [ex, ey, ez])
        return 0., 0., 0.
    end

    dlat_dx   = -z / (x^2 + y^2)^(3 / 2) / (1 + z^2 / (x^2 + y^2)) * x
    dlat_dy   = -z / (x^2 + y^2)^(3 / 2) / (1 + z^2 / (x^2 + y^2)) * y
    dlat_dz   = 1. / (x^2 + y^2)^(1 / 2) / (1 + z^2 / (x^2 + y^2))
    dlon_dx   = -y / x^2 / (1. + (y / x)^2)
    dlon_dy   = 1 / x / (1. + (y / x)^2)
    dlon_dz   = 0.
    dr_dx   = x / sqrt(x^2 + y^2 + z^2)
    dr_dy   = y / sqrt(x^2 + y^2 + z^2)
    dr_dz   = z / sqrt(x^2 + y^2 + z^2)

    M =  [dlon_dx dlon_dy dlon_dz;
          dlat_dx dlat_dy dlat_dz;
          dr_dx   dr_dy   dr_dz]

    elon, elat, erotrate = sqrt.(diag(M * diagm([ex, ey, ez].^2) * M'))
    erotrate = 1e6 * rad2deg(erotrate)

    elon, elat, erotrate
end


function err_cart_to_sphere(pole::PoleCart)
    err_cart_to_sphere(pole.x, pole.y, pole.z, pole.ex, pole.ey, pole.ez)
end


function err_cart_to_sphere(x::Float64, y::Float64, z::Float64, 
                            ex::Float64, ey::Float64, ez::Float64)

    if all(x -> x == 0., [ex, ey, ez])
        return 0., 0., 0.
    end
    
    r = sqrt(x^2 + y^2 + z^2)
    
    x_norm = x / r
    y_norm = y / r
    z_norm = z / r

    theta = atan(y, x)
    phi = atan(z_norm / sqrt(x_norm^2 + y_norm^2))

    d11 = cos(theta) * sin(phi)
    d12 = sin(theta) * sin(phi)
    d13 = cos(phi)

    d21 = - sin(theta) / (r * sin(phi))
    d22 = cos(theta) / (r * sin(phi))
    d23 = 0.
    
    d31 = cos(theta) * cos(phi) / r
    d32 = sin(theta) * cos(phi) / r
    d33 = - sin(phi) / r

    M = [d11 d12 d13; d21 d22 d23; d31 d32 d33]

    er, etheta, ephi = abs.(M * [ex, ey, ez])
    # er, etheta, ephi = sqrt.(diag(M * diagm([ex, ey, ez].^2) * M'))

    elon = rad2deg(etheta)
    elat = rad2deg(ephi)
    erotrate = 1e6 * rad2deg(er)

    elon, elat, erotrate
end


function err_sphere_to_cart(lon::Float64, lat::Float64, rotrate::Float64, 
                            elon::Float64, elat::Float64, erotrate::Float64)

    if all(x -> x == 0., [elon, elat, erotrate])
        return 0., 0., 0.
    end

    theta = deg2rad(lon)
    phi = deg2rad(lat)
    r = deg2rad(rotrate) / 1e6

    etheta = deg2rad(elon)
    ephi = deg2rad(elat)
    er = deg2rad(erotrate) / 1e6

    d11 = cos(theta) * sin(phi)
    d12 = -r * sin(theta) * sin(phi)
    d13 = r * cos(theta) * cos(phi)

    d21 = sin(theta) * sin(phi)
    d22 = r * cos(theta) * sin(phi)
    d23 = r * sin(theta) * cos(phi)

    d31 = cos(phi)
    d32 = 0.
    d33 = - r * sin(phi)

    M = [d11 d12 d13; d21 d22 d23; d31 d32 d33]

    # ex, ey, ez = sqrt.(diag(M * diagm([er, etheta, ephi].^2) * M'))
    ex, ey, ez = abs.(M * [er, etheta, ephi])

end


function err_sphere_to_cart(pole::PoleSphere)
    err_sphere_to_cart(pole.lon, pole.lat, pole.rotrate,
                       pole.elon, pole.elat, pole.erotrate)
end




"""
    pole_cart_to_sphere(pole)

Converts an `PoleCart` (an Euler vector in Cartesian coordinates)
into spherical coordinates (lon, lat, deg / Myr).
"""
function pole_cart_to_sphere(pole::PoleCart)
    rotation_rate_cart = sqrt(pole.x^2 + pole.y^2 + pole.z^2)
    
    pole_x_norm = pole.x / rotation_rate_cart
    pole_y_norm = pole.y / rotation_rate_cart
    pole_z_norm = pole.z / rotation_rate_cart

    pole_lon = atand(pole_y_norm, pole_x_norm)
    pole_lat = atand(pole_z_norm / sqrt(pole_x_norm^2 + pole_y_norm^2))

    rotation_rate_deg_Myr = rad2deg(rotation_rate_cart) * 1e6

    elon, elat, erotrate = err_cart_to_sphere(pole)

    PoleSphere(lon=pole_lon, lat=pole_lat, rotrate=rotation_rate_deg_Myr,
               elon=elon, elat=elat, erotrate=erotrate,
               fix=pole.fix, mov=pole.mov)
end


function pole_sphere_to_cart(pole::PoleSphere)
    r = deg2rad(pole.rotrate) / 1e6

    x = r * cosd(pole.lat) * cosd(pole.lon)
    y = r * cosd(pole.lat) * sind(pole.lon)
    z = r * sind(pole.lat)

    ex, ey, ez = err_sphere_to_cart(pole)

    PoleCart(x=x, y=y, z=z, ex=ex, ey=ey, ez=ez, fix=pole.fix, mov=pole.mov)
end


end

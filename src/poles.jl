module RotationPoles

export PoleCart, PoleSphere, add_poles, add_poles, subtract_poles,
    pole_cart_to_sphere, pole_sphere_to_cart

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
    fix::String = ""
    mov::String = ""
end


function Base.:-(pole::PoleCart)
    PoleCart(x = -pole.x, y = -pole.y, z = -pole.z, fix = pole.mov, mov = pole.fix)
end

function Base.:-(pole::PoleSphere)
    pc = pole_sphere_to_cart(pole)
    pole_cart_to_sphere(-pc)
end


function add_poles(poles::PoleCart...)
    xn = sum(pole.x for pole in poles)
    yn = sum(pole.y for pole in poles)
    zn = sum(pole.z for pole in poles)

    PoleCart(x = xn, y = yn, z = zn, fix = poles[1].fix, 
                  mov = poles[end].mov)
end


function add_poles(poles::Array{PoleCart})
    xn = sum(pole.x for pole in poles)
    yn = sum(pole.y for pole in poles)
    zn = sum(pole.z for pole in poles)

    PoleCart(x = xn, y = yn, z = zn, fix = poles[1].fix, 
                  mov = poles[end].mov)
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

    PoleCart(x = xx, y = yy, z = zz, fix = pole1.fix, mov = pole2.fix)
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

"""
    euler_pole_cart_to_sphere(pole)

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

    PoleSphere(lon = pole_lon, lat = pole_lat, 
                    rotrate = rotation_rate_deg_Myr,
    fix = pole.fix, mov = pole.mov)
end


function pole_sphere_to_cart(pole::PoleSphere)
    r = deg2rad(pole.rotrate) / 1e6

    x = r * cosd(pole.lat) * cosd(pole.lon)
    y = r * cosd(pole.lat) * sind(pole.lon)
    z = r * sind(pole.lat)

    PoleCart(x = x, y = y, z = z, fix = pole.fix, mov = pole.mov)
end


end

"""
    VelocityVectorSph(20., )
Velocity vector in spherical (lon, lat) coordinates, with velocities in
mm/yr.

ve, vn, vu are east, north and up velocities, and ee, en, and eu are the 
1-sigma uncertainties.
"""
struct VelocityVectorSph
    latd::AbstractFloat
    lond::AbstractFloat
    ve::AbstractFloat
    vn::AbstractFloat
    vu::AbstractFloat
    ee::AbstractFloat
    en::AbstractFloat
    eu::AbstractFloat
end


function build_Pv_deg(lond::AbstractFloat, latd::AbstractFloat)

    pv11 = -sind(latd) * cosd(lond)
    pv12 = -sind(latd) * sind(lond)
    pv13 = cosd(latd)

    pv21 = -sind(lond)
    pv22 = cosd(lond)
    pv23 = 0

    pv31 = -cosd(latd) * cosd(lond)
    pv32 = -cosd(latd) * sind(lond)
    pv33 = -sind(latd)

    Pv = [pv11 pv12 pv13;
          pv12 pv22 pv23;
          pv31 pv32 pv33]
    return Pv
end


function build_Gb_deg(lond::AbstractFloat, latd::AbstractFloat; R=6371000.)
    x_hat = R * cosd(latd) * cosd(lond)
    y_hat = R * cosd(latd) * sind(lond)
    z_hat = R * sind(latd);

    Gb = [ 0      z_hat  -y_hat;
          -z_hat  0       x_hat;
           y_hat -x_hat   0]
end


"""
    build_PvGb_deg(lond, latd)

Creates a linear operator to convert a Cartesian Euler vector into
a spherical velocity vector at a given (lon, lat) in degrees; in other
words, this predicts the velocity of a point on the earth's surface
(such as that measured by a GPS station or from a fault) given an Euler
pole.

# Arguments
- `lond::AbstractFloat`: The longitude of the point, in degrees
- `latd::AbstractFloat`: The latitude of the point, in degrees

# Returns
- `PvGb`: 3x3 AbstractFloat matrix

# Examples
```jldoctest
julia> pg = build_PvGb_deg(0., 0.)
3x3 Array{Float64,2}:
 0.0  -6.371e6  0.0
 0.0   0.0      6.371e6
 0.0   0.0      0.0
```
"""
function build_PvGb_deg(lond::AbstractFloat, latd::AbstractFloat)
    Pv = build_Pv_deg(lond, latd)
    Gb = build_Gb_deg(lond, latd)

    PvGb = Pv * Gb
end


"""
    build_PvGb_vel(vel)

Creates a linear operator to convert a Cartesian Euler vector into
a spherical velocity vector at a given (lon, lat) in degrees; in other
words, this predicts the velocity of a point on the earth's surface
(such as that measured by a GPS station or from a fault) given an Euler
pole.

# Arguments
- `vel::VelocityVectorSph`: A velocity vector with `lond` and `latd` attributes

# Returns
- `PvGb`: 3x3 AbstractFloat matrix

# Examples
```jldoctest
julia> pg = build_PvGb_deg(0., 0.)
3x3 Array{Float64,2}:
 0.0  -6.371e6  0.0
 0.0   0.0      6.371e6
 0.0   0.0      0.0
```
"""
function build_PvGb_vel(vel::VelocityVectorSph)
    build_PvGb_deg(vel.lond, vel.latd)
end


"""
"""
function build_PvGb_from_vels(vels::Array{VelocityVectorSph,1})
    reduce(vcat, [build_PvGb_vel(vel) for vel in vels])
end


function build_vel_column_from_vel(vel::VelocityVectorSph)
    V = [vel.ve; vel.vn; vel.vu]
end


function build_vel_column_from_vels(vels::Array{VelocityVectorSph,1})
    reduce(vcat, [build_vel_column_from_vel(vel) for vel in vels])
end


"""
Euler Pole (rotation vector) in Cartesian coordinates
"""
struct EulerPoleCart
    x::AbstractFloat
    y::AbstractFloat
    z::AbstractFloat
end


"""
Euler Pole (rotation vector) in spherical coordinates.
"""
struct EulerPoleSphere
    lond::AbstractFloat
    latd::AbstractFloat
    rotrate::AbstractFloat
end


function add_poles(pole1::EulerPoleCart, pole2::EulerPoleCart)
    xn = pole1.x + pole2.x
    yn = pole1.y + pole2.y
    zn = pole1.z + pole2.z

    EulerPoleCart(xn, yn, zn)
end


function add_poles(pole1::EulerPoleSphere, pole2::EulerPoleSphere)
    cart_pole_1 = euler_pole_sphere_to_cart(pole1)
    cart_pole_2 = euler_pole_sphere_to_cart(pole2)

    cart_pole_sum = add_poles(cart_pole_1, cart_pole_2)
    euler_pole_cart_to_sphere(cart_pole_sum)
end


"""
    euler_pole_cart_to_sphere(pole)

Converts an `EulerPoleCart` (an Euler vector in Cartesian coordinates)
into spherical coordinates (lond, latd, deg / Myr).
"""
function euler_pole_cart_to_sphere(pole::EulerPoleCart)
    rotation_rate_cart = sqrt(pole.x^2 + pole.y^2 + pole.z^2)
    
    pole_x_norm = pole.x / rotation_rate_cart
    pole_y_norm = pole.y / rotation_rate_cart
    pole_z_norm = pole.z / rotation_rate_cart

    pole_lon = atand(pole_y_norm, pole_x_norm)
    pole_lat = acosd(pole_z_norm)

    rotation_rate_deg_Myr = rad2deg(rotation_rate_cart) * 1e6

    EulerPoleSphere(pole_lon, pole_lat, rotation_rate_deg_Myr)
end


function euler_pole_sphere_to_cart(pole::EulerPoleSphere)
    r = deg2rad(pole.rotrate / 1e6)

    x = r * sind(pole.latd) * cosd(pole.lond)
    y = r * sind(pole.latd) * sind(pole.lond)
    z = r * cosd(pole.latd)

    EulerPoleCart(x, y, z)
end
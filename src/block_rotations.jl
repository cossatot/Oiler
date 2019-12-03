using Parameters

"""
    VelocityVectorSph(20., )
Velocity vector in spherical (lon, lat) coordinates, with velocities in
mm/yr.

ve, vn, vu are east, north and up velocities, and ee, en, and eu are the 
1-sigma uncertainties.
"""
@with_kw struct VelocityVectorSph
    lond::Float64
    latd::Float64
    ve::Float64
    vn::Float64
    vu::Float64 = 0.
    ee::Float64 = 0.
    en::Float64 = 0.
    eu::Float64 = 0.
end


function build_Pv_deg(lond::Float64, latd::Float64)

    pv11 = -sind(latd) * cosd(lond)
    pv12 = -sind(latd) * sind(lond)
    pv13 = cosd(latd)

    pv21 = -sind(lond)
    pv22 = cosd(lond)
    pv23 = 0.

    pv31 = -cosd(latd) * cosd(lond)
    pv32 = -cosd(latd) * sind(lond)
    pv33 = -sind(latd)

    Pv = [pv11 pv12 pv13;
          pv21 pv22 pv23;
          pv31 pv32 pv33]
    return Pv
end


function build_Gb_deg(lond::Float64, latd::Float64; R = 6371000.)
    x_hat = R * cosd(latd) * cosd(lond)
    y_hat = R * cosd(latd) * sind(lond)
    z_hat = R * sind(latd);

    Gb = [ 0.     z_hat  -y_hat;
          -z_hat  0.      x_hat;
           y_hat -x_hat   0.]
end


"""
    build_PvGb_deg(lond, latd)

Creates a linear operator to convert a Cartesian Euler vector into
a spherical velocity vector at a given (lon, lat) in degrees; in other
words, this predicts the velocity of a point on the earth's surface
(such as that measured by a GPS station or from a fault) given an Euler
pole.

# Arguments
- `lond::Float64`: The longitude of the point, in degrees
- `latd::Float64`: The latitude of the point, in degrees

# Returns
- `PvGb`: 3x3 Float64 matrix

# Examples
```jldoctest
julia> pg = build_PvGb_deg(0., 0.)
3x3 Array{Float64,2}:
 0.0  -6.371e6  0.0
 0.0   0.0      6.371e6
 0.0   0.0      0.0
```
"""
function build_PvGb_deg(lond::Float64, latd::Float64)
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
- `PvGb`: 3x3 Float64 matrix

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
function build_PvGb_from_vels(vels::Array{VelocityVectorSph})
    reduce(vcat, [build_PvGb_vel(vel) for vel in vels])
end


function build_PvGb_from_degs(londs::Array{Float64},
                              latds::Array{Float64})
    reduce(vcat, [build_PvGb_deg(lond, latds[i])
                  for (i, lond) in enumerate(londs)])
end

function build_vel_column_from_vel(vel::VelocityVectorSph)
    V = [vel.ve; vel.vn; vel.vu]
end


function build_vel_column_from_vels(vels::Array{VelocityVectorSph})
    reduce(vcat, [build_vel_column_from_vel(vel) for vel in vels])
end


"""
Euler Pole (rotation vector) in Cartesian coordinates

# Arguments
- `x::Float64`: *x* coordinate in the Cartesian system
- `y::Float64`: *y* coordinate in the Cartesian system
- `z::Float64`: *z* coordinate in the Cartesian system
- `fix::String`: Name of fixed block. Defaults to "".
- `rel::String`: Name of moving (relative) block.  Defaults to "".
"""
@with_kw struct EulerPoleCart
    x::Float64
    y::Float64
    z::Float64
    fix::String = ""
    rel::String = ""
end


"""
Euler Pole (rotation vector) in spherical coordinates.

# Arguments
- `lond::Float64`: Longitude coordinate in the spherical system (in degrees)
- `latd::Float64`: Latitude coordinate in the spherical system (in degrees)
- `rotrate::Float64`: Rotation rate (in degree / time unit, unspecified)
- `fix::String`: Name of fixed block. Defaults to "".
- `rel::String`: Name of moving (relative) block.  Defaults to "".
"""
@with_kw struct EulerPoleSphere
    lond::Float64
    latd::Float64
    rotrate::Float64
    fix::String = ""
    rel::String = ""
end


function add_poles(poles::EulerPoleCart...)
    xn = sum(pole.x for pole in poles)
    yn = sum(pole.y for pole in poles)
    zn = sum(pole.z for pole in poles)

    EulerPoleCart(x = xn, y = yn, z = zn, fix = poles[1].fix, 
                  rel = poles[end].rel)
end


function add_poles(poles::EulerPoleSphere...)
    cart_poles = [euler_pole_sphere_to_cart(pole) for pole in poles]

    cart_pole_sum = add_poles(cart_poles...)
    euler_pole_cart_to_sphere(cart_pole_sum)
end


function subtract_poles(pole1::EulerPoleCart, pole2::EulerPoleCart)
    xx = pole1.x - pole2.x
    yy = pole1.y - pole2.y
    zz = pole1.z - pole2.z

    EulerPoleCart(x = xx, y = yy, z = zz, fix = pole1.fix, rel = pole2.fix)
end

function subtract_poles(pole1::EulerPoleSphere, pole2::EulerPoleSphere)
    pole1c = euler_pole_sphere_to_cart(pole1)
    pole2c = euler_pole_sphere_to_cart(pole2)

    pole_diff = subtract_poles(pole1c, pole2c)

    euler_pole_cart_to_sphere(pole_diff)
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
    pole_lat = atand(pole_z_norm / sqrt(pole_x_norm^2 + pole_y_norm^2))

    rotation_rate_deg_Myr = rad2deg(rotation_rate_cart) * 1e6

    EulerPoleSphere(lond = pole_lon, latd = pole_lat, 
                    rotrate = rotation_rate_deg_Myr,
    fix = pole.fix, rel = pole.rel)
end


function euler_pole_sphere_to_cart(pole::EulerPoleSphere)
    r = deg2rad(pole.rotrate / 1e6)

    x = r * cosd(pole.latd) * cosd(pole.lond)
    y = r * cosd(pole.latd) * sind(pole.lond)
    z = r * sind(pole.latd)

    EulerPoleCart(x = x, y = y, z = z, fix = pole.fix, rel = pole.rel)
end


function predict_block_vels(londs::Array{Float64},
                            latds::Array{Float64},
                            pole::EulerPoleSphere)

    PvGb = build_PvGb_from_degs(londs, latds)

    cart_pole = euler_pole_sphere_to_cart(pole)

    V_pred = PvGb * [cart_pole.x; cart_pole.y; cart_pole.z]
    Ve_pred = V_pred[1:3:end]
    Vn_pred = V_pred[2:3:end]
    Vu_pred = V_pred[3:3:end]

    n_vels = length(londs)

    pred_vels = Array{VelocityVectorSph}(undef, n_vels)

    for n in 1:n_vels
        pred_vels[n] = VelocityVectorSph(londs[n], latds[n],
                                         Ve_pred[n], Vn_pred[n], Vu_pred[n],
                                         0., 0., 0.)
    end
    return pred_vels
end
    
function predict_block_vels(londs::Array{Float64},
                            latds::Array{Float64},
                            pole::EulerPoleCart)

    PvGb = build_PvGb_from_degs(londs, latds)

    V_pred = PvGb * [pole.x; pole.y; pole.z]
    Ve_pred = V_pred[1:3:end]
    Vn_pred = V_pred[2:3:end]
    Vu_pred = V_pred[3:3:end]

    if length(size(londs)) == 1
        n_vels = size(londs)[1]
    else
        n_vels = size(londs)[2]
    end

    pred_vels = Array{VelocityVectorSph}(undef, n_vels)

    for n in 1:n_vels
        pred_vels[n] = VelocityVectorSph(londs[n], latds[n],
                                         Ve_pred[n], Vn_pred[n], Vu_pred[n],
                                         0., 0., 0.)
    end
    return pred_vels
end
    


function calc_strike(lond1::Float64, latd1::Float64,
                     lond2::Float64, latd2::Float64)
        
    y = sind(lond2 - lond1) * cosd(latd1)
    x = cosd(latd1) * sind(latd2) - sind(latd1) * cosd(latd2) * cos(lond2 - lond1)

    strike = atand(y, x)
end
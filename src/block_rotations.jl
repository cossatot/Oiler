module BlockRotations

export build_PvGb_from_vels, build_vel_column_from_vels, predict_block_vels

using LinearAlgebra
using ..Oiler: EARTH_RAD_MM, VelocityVectorSphere, PoleSphere, PoleCart,
        pole_sphere_to_cart, pole_cart_to_sphere


function build_Pv_deg(lon::Float64, lat::Float64)
    pex = -sind(lon)
    pey = cosd(lon)
    pez = 0

    pnx = -sind(lat) * cosd(lon)
    pny = -sind(lat) * sind(lon)
    pnz = cosd(lat)

    pux = cosd(lat) * cosd(lon)
    puy = cosd(lat) * sind(lon)
    puz = sind(lat)

    Pv = [pex pey pez;
          pnx pny pnz;
          pux puy puz]
end


function build_Gb_deg(lon::Float64, lat::Float64; R=EARTH_RAD_MM)
    x_hat = R * cosd(lat) * cosd(lon)
    y_hat = R * cosd(lat) * sind(lon)
    z_hat = R * sind(lat);

    Gb = [ 0.     z_hat  -y_hat;
          -z_hat  0.      x_hat;
           y_hat -x_hat   0.]
end


"""
    build_PvGb_deg(lon, lat)

Creates a linear operator to convert a Cartesian Euler vector into
a spherical velocity vector at a given (lon, lat) in degrees; in other
words, this predicts the velocity of a point on the earth's surface
(such as that measured by a GPS station or from a fault) given an Euler
pole.

# Arguments
- `lon::Float64`: The longitude of the point, in degrees
- `lat::Float64`: The latitude of the point, in degrees

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
function build_PvGb_deg(lon::Float64, lat::Float64)
    Pv = build_Pv_deg(lon, lat)
    Gb = build_Gb_deg(lon, lat)

    PvGb = Pv * Gb
end


"""
    build_PvGb_vel(vel::VelocityVectorSphere)

Creates a linear operator to convert a Cartesian Euler vector into
a spherical velocity vector at a given (lon, lat) in degrees; in other
words, this predicts the velocity of a point on the earth's surface
(such as that measured by a GPS station or from a fault) given an Euler
pole.

# Arguments
- `vel::VelocityVectorSphere`: A velocity vector with `lon` and `lat` attributes

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
function build_PvGb_vel(vel::VelocityVectorSphere)
    build_PvGb_deg(vel.lon, vel.lat)
end


"""

"""
function build_PvGb_from_vels(vels::Array{VelocityVectorSphere})
    reduce(vcat, [build_PvGb_vel(vel) for vel in vels])
end


function build_PvGb_from_degs(lons::Array{Float64},
                              lats::Array{Float64})
    reduce(vcat, [build_PvGb_deg(lon, lats[i])
                  for (i, lon) in enumerate(lons)])
end

function build_vel_column_from_vel(vel::VelocityVectorSphere)
    V = [vel.ve; vel.vn; vel.vu]
end


function build_vel_column_from_vels(vels::Array{VelocityVectorSphere})
    reduce(vcat, [build_vel_column_from_vel(vel) for vel in vels])
end


function predict_block_vel(lon::Float64, lat::Float64, pole::PoleCart)
    vel_pred = predict_block_vels([lon], [lat], pole)[1]
end
    

function predict_block_vel(lon::Float64, lat::Float64, pole::PoleSphere)
    predict_block_vel(lon, lat, pole_sphere_to_cart(pole))
end


function predict_block_vels(lons::Array{Float64},
                            lats::Array{Float64},
                            pole::PoleSphere)

    PvGb = build_PvGb_from_degs(lons, lats)

    cart_pole = pole_sphere_to_cart(pole)

    V_pred = PvGb * [cart_pole.x; cart_pole.y; cart_pole.z]
    Ve_pred = V_pred[1:3:end]
    Vn_pred = V_pred[2:3:end]
    Vu_pred = V_pred[3:3:end]

    n_vels = length(lons)

    pred_vels = Array{VelocityVectorSphere}(undef, n_vels)

    for n in 1:n_vels
        pred_vels[n] = VelocityVectorSphere(lon=lons[n], lat=lats[n],
                                            ve=Ve_pred[n], vn=Vn_pred[n],
                                            fix=pole.fix, mov=pole.mov)
    end
    return pred_vels
end
    

function predict_block_vels(lons::Array{Float64},
                            lats::Array{Float64},
                            pole::PoleCart)

    PvGb = build_PvGb_from_degs(lons, lats)

    V_pred = PvGb * [pole.x; pole.y; pole.z]
    Ve_pred = V_pred[1:3:end]
    Vn_pred = V_pred[2:3:end]
    Vu_pred = V_pred[3:3:end]

    if length(size(lons)) == 1
        n_vels = size(lons)[1]
    else
        n_vels = size(lons)[2]
    end

    # Propagate uncertainty in pole (Vel locations have negligible uncertainty)
    if any(x -> x != 0., [pole.ex, pole.ey, pole.ez])
        V_err = diag(PvGb * diagm([pole.ex, pole.ey, pole.ez]) * PvGb')
    else
        V_err = zeros(size(V_pred))
    end
    Ve_err = V_err[1:3:end]
    Vn_err = V_err[2:3:end]

    pred_vels = Array{VelocityVectorSphere}(undef, n_vels)

    for n in 1:n_vels
        pred_vels[n] = VelocityVectorSphere(lon=lons[n], lat=lats[n],
                                            ve=Ve_pred[n], vn=Vn_pred[n],
                                            ee=Ve_err[n], en=Vn_err[n],
                                            fix=pole.fix, mov=pole.mov)
    end
    return pred_vels
end
   

function predict_block_vels(vels::Array{VelocityVectorSphere},
pole::PoleCart)
    lons = [v.lon for v in vels]
    lats = [v.lat for v in vels]

    predict_block_vels(lons, lats, pole)
end


function predict_block_vels(vels::Array{VelocityVectorSphere}, 
    pole::PoleSphere)
    predict_block_vels(vels, pole_sphere_to_cart(pole))
end

end # module
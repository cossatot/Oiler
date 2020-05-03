module BlockRotations

export build_PvGb_from_vels, build_vel_column_from_vels, predict_block_vels

using ..Oiler: EARTH_RAD_MM, VelocityVectorSphere, PoleSphere, PoleCart,
        pole_sphere_to_cart, pole_cart_to_sphere


function build_Pv_deg_ned(lond::Float64, latd::Float64)

    pv11 = -sind(latd) * cosd(lond) #nx
    pv12 = -sind(latd) * sind(lond) #ny
    pv13 = cosd(latd) #nz

    pv21 = -sind(lond) #ex
    pv22 = cosd(lond)  #ey
    pv23 = 0. #ez

    pv31 = -cosd(latd) * cosd(lond) #dx
    pv32 = -cosd(latd) * sind(lond) #dy
    pv33 = -sind(latd) #dz

    Pv = [pv11 pv12 pv13;
          pv21 pv22 pv23;
          pv31 pv32 pv33]
    return Pv
end


function build_Pv_deg(lond::Float64, latd::Float64)
    pex = -sind(lond)
    pey = cosd(lond)
    pez = 0

    pnx = -sind(latd) * cosd(lond)
    pny = -sind(latd) * sind(lond)
    pnz = cosd(latd)

    pux = cosd(latd) * cosd(lond)
    puy = cosd(latd) * sind(lond)
    puz = sind(latd)

    Pv = [pex pey pez;
          pnx pny pnz;
          pux puy puz]
end


function build_Gb_deg(lond::Float64, latd::Float64; R = EARTH_RAD_MM)
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
    build_PvGb_vel(vel::VelocityVectorSphere)

Creates a linear operator to convert a Cartesian Euler vector into
a spherical velocity vector at a given (lon, lat) in degrees; in other
words, this predicts the velocity of a point on the earth's surface
(such as that measured by a GPS station or from a fault) given an Euler
pole.

# Arguments
- `vel::VelocityVectorSphere`: A velocity vector with `lond` and `latd` attributes

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
    build_PvGb_deg(vel.lond, vel.latd)
end


"""

"""
function build_PvGb_from_vels(vels::Array{VelocityVectorSphere})
    reduce(vcat, [build_PvGb_vel(vel) for vel in vels])
end


function build_PvGb_from_degs(londs::Array{Float64},
                              latds::Array{Float64})
    reduce(vcat, [build_PvGb_deg(lond, latds[i])
                  for (i, lond) in enumerate(londs)])
end

function build_vel_column_from_vel(vel::VelocityVectorSphere)
    V = [vel.ve; vel.vn; vel.vu]
end


function build_vel_column_from_vels(vels::Array{VelocityVectorSphere})
    reduce(vcat, [build_vel_column_from_vel(vel) for vel in vels])
end



function predict_block_vels(londs::Array{Float64},
                            latds::Array{Float64},
                            pole::PoleSphere)

    PvGb = build_PvGb_from_degs(londs, latds)

    cart_pole = pole_sphere_to_cart(pole)

    V_pred = PvGb * [cart_pole.x; cart_pole.y; cart_pole.z]
    Ve_pred = V_pred[1:3:end]
    Vn_pred = V_pred[2:3:end]
    Vu_pred = V_pred[3:3:end]

    n_vels = length(londs)

    pred_vels = Array{VelocityVectorSphere}(undef, n_vels)

    for n in 1:n_vels
        pred_vels[n] = VelocityVectorSphere(lond = londs[n], latd = latds[n],
                                            ve = Ve_pred[n], vn = Vn_pred[n],
                                            fix = pole.fix, mov = pole.mov)
    end
    return pred_vels
end
    

function predict_block_vels(londs::Array{Float64},
                            latds::Array{Float64},
                            pole::PoleCart)

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

    pred_vels = Array{VelocityVectorSphere}(undef, n_vels)

    for n in 1:n_vels
        pred_vels[n] = VelocityVectorSphere(lond = londs[n], latd = latds[n],
                                            ve = Ve_pred[n], vn = Vn_pred[n],
                                            fix = pole.fix, mov = pole.mov)
    end
    return pred_vels
end
   

function predict_block_vels(vels::Array{VelocityVectorSphere},
pole::PoleCart)
    londs = [v.lond for v in vels]
    latds = [v.latd for v in vels]

    predict_block_vels(londs, latds, pole)
end


function predict_block_vels(vels::Array{VelocityVectorSphere}, 
    pole::PoleSphere)
    predict_block_vels(vels, pole_sphere_to_cart(pole))
end

end #module
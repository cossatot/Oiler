module Faults
export Fault, fault_to_vel, fault_slip_rate_to_ve_vn, ve_vn_to_fault_slip_rate,
    fault_oblique_merc

using Parameters

using ..Oiler: VelocityVectorSphere, average_azimuth, az_to_angle,
angle_difference, rotate_velocity, EARTH_RAD_KM

const direction_map = Dict{String,Float64}("N" => 0.,
                 "NNE" => 22.5,
                 "NE" => 45.,
                 "ENE" => 67.5,
                 "E" => 90.,
                 "ESE" => 112.5,
                 "SE" => 135.,
                 "SSE" => 157.5,
                 "S" => 180.,
                 "SSW" => 202.5,
                 "SW" => 225.,
                 "WSW" => 247.5,
                 "W" => 270.,
                 "WNW" => 292.5,
                 "NW" => 315.,
                 "NNW" => 337.5)


"""
    Fault(trace, strike, dip, dip_dir, extension_rate, extension_err,
        dextral_rate, dextral_err, lsd, usd, name, hw, fw)
Representation of the geometric and kinematic attributes of a fault.

# Arguments:
 
- `trace::Array{Float64,2}`: Surface coordinates of the trace. Rows are each
    point, with the longitude in the first column and latitude in the second
    column.

- `strike::Float64`: Strike (degrees, azimuth) of the trace in right-hand-rule 
    system.

- `dip::Float64`: Dip of the fault in degrees.

- `dip_dir::String`: Cardinal direction that the fault dips towards. Acceptable
    values are in capital letters and may represent primary, secondary or
    tertiary cardinal direction, e.g. "N", "NE", "ENE". Vertically-dipping
    faults should have a `dip_dir` of "V".

- `extension_rate::Float64 = 0.`: Rate of horizontal extension across the fault
    in mm/yr.  Negative values indicate shortening.

- `extension_err:Float64 = 0.`: 1-s.d. uncertainty of horizontal extension rate
    in mm/yr.

- `dextral_rate::Float64 = 0.`: Rate of dextral slip on the fault in mm/yr.
    Negative values indicate sinistral slip.
- `dextral_err::Float64 = 0.`: 1-s.d. uncertainty of dextral slip rate in mm/yr.
- `lsd::Float64 = 25.`: Lower seismogenic/locking depth, in km.
- `usd::Float64 = 0.`: Upper seismogenic/locking_depth, in km.
- `name::String = ""`: Name of the fault (segment).
- `hw::String = ""`: Name of the hanging wall block; for vertically-dipping
    faults, ...
- `fw::String = ""`: Name of the footwall block; for vertically-dipping faults,
    ...
"""
struct Fault
    trace::Array{Float64,2}
    strike::Float64
    dip::Float64
    dip_dir::String
    extension_rate::Float64
    extension_err::Float64
    dextral_rate::Float64
    dextral_err::Float64
    lsd::Float64
    usd::Float64
    name::String
    hw::String
    fw::String
end


"""
    fault(trace, dip, dip_dir, extension_rate, extension_err,
        dextral_rate, dextral_err, lsd, usd, name, hw, fw)
Constructor function for the `Fault` struct that represents the geometric and 
kinematic attributes of a fault.  The strike is calculated from the trace
coordinates and the dip direction to ensure right-hand-rule compatibility.

# Arguments:
 
- `trace::Array{Float64,2}`: Surface coordinates of the trace. Rows are each
    point, with the longitude in the first column and latitude in the second
    column.

- `dip::Float64`: Dip of the fault in degrees.

- `dip_dir::String`: Cardinal direction that the fault dips towards. Acceptable
    values are in capital letters and may represent primary, secondary or
    tertiary cardinal direction, e.g. "N", "NE", "ENE".

- `extension_rate::Float64 = 0.`: Rate of horizontal extension across the fault
    in mm/yr.  Negative values indicate shortening.

- `extension_err:Float64 = 0.`: 1-s.d. uncertainty of horizontal extension rate
    in mm/yr.

- `dextral_rate::Float64 = 0.`: Rate of dextral slip on the fault in mm/yr.
    Negative values indicate sinistral slip.
- `dextral_err::Float64 = 0.`: 1-s.d. uncertainty of dextral slip rate in mm/yr.
- `lsd::Float64 = 25.`: Lower seismogenic/locking depth, in km.
- `usd::Float64 = 0.`: Upper seismogenic/locking_depth, in km.
- `name::String = ""`: Name of the fault (segment).
- `hw::String = ""`: Name of the hanging wall block; for vertically-dipping
    faults, this is the block to the left-hand side when looking from the last
    coordinate to the first coordinate.
- `fw::String = ""`: Name of the footwall block; for vertically-dipping faults,
    this is the block to the right-hand side when looking from the last
    coordinate to the first coordinate.
    ...
"""
function Fault(; trace::Array{Float64,2}, dip::Float64,
    dip_dir::String,
    extension_rate::Float64 = 0.,
    extension_err::Float64 = 0.,
    dextral_rate::Float64 = 0.,
    dextral_err::Float64 = 0.,
    lsd::Float64 = 25.,
    usd::Float64 = 0.,
    name::String =  "",
    hw::String = "",
    fw::String = "")

    if dip != 90.
        trace = check_right_hand_rule(trace, dip_dir)
    end

    strike = average_azimuth(trace[:,1], trace[:,2])

    # Fault(trace = trace, strike = strike, dip = dip, dip_dir = dip_dir, 
    #    extension_rate = extension_rate, extension_err = extension_err, 
    #    dextral_rate = dextral_rate, dextral_err = dextral_err, 
    #    lsd = lsd, usd = usd, name = name, hw = hw, fw = fw)
    Fault(trace, strike, dip, dip_dir, 
        extension_rate, extension_err, 
        dextral_rate, dextral_err, 
        lsd, usd, name, hw, fw)
end

"""
    fault_to_vel(fault)
Converts a `Fault` object into a `VelocityVectorSphere`.  

The velocity is expressed as the horizontal velocity of the moving footwall wall
relative to the fixed hanging wall (for vertical faults, this is the moving
right-hand-side block relative to the fixed left-hand-side block when looking
from the last trace coordinate to the first trace coordinate). The geographic
location of the `VelocityVectorSphere` is either the middle trace coordinate (if
there is an odd number of vertices) or the middle of the central section (if
there is an even number of vertices); this may or may not correspond to the
exact midpoint of the fault trace, depending on the distribution of vertices
along the fault's length.
"""
function fault_to_vel(fault::Fault)
    ve, vn = fault_slip_rate_to_ve_vn(fault.dextral_rate, fault.extension_rate,
        fault.strike)
    
    # this may not be valid, need to check...
    ee, en = fault_slip_rate_to_ve_vn(fault.dextral_err, fault.extension_err,
        fault.strike)

    vlon, vlat = get_midpoint(fault.trace)

    VelocityVectorSphere(lond = vlon, latd = vlat, ve = ve, vn = vn, 
        fix = fault.hw, mov = fault.fw, name = fault.name, ee = ee, en = en,
        vel_type = "fault")
end


function get_midpoint(trace::Array{Float64,2})
    n_pts = size(trace)[1]

    if iseven(n_pts)
        mid = n_pts ÷ 2
        return ((trace[mid, 1] + trace[mid + 1, 1]) / 2.,
                (trace[mid, 2] + trace[mid + 1, 2]) / 2.)
    else
        mid = n_pts ÷ 2 + 1
        return (trace[mid, 1], trace[mid, 2])
    end
end


function check_right_hand_rule(trace::Array{Float64,2}, dip_dir::String; 
    reverse_angle_threshold::Float64 = 90.)
    # Modified from the OQ-MBTK tools, (c) Global Earthquake Model Foundation

    strike = average_azimuth(trace[:,1], trace[:,2])

    trace_dip_trend = strike + 90.

    fault_dip_trend = direction_map[dip_dir]

    trend_angle_difference = angle_difference(trace_dip_trend, fault_dip_trend)

    # TODO: warn if trend_angle_difference < 15 (or some other low threshold)
    if trend_angle_difference > reverse_angle_threshold
        trace = reverse(trace, dims = 1)
    end
    trace
end


function fault_slip_rate_to_ve_vn(dextral_rate::Float64, extension_rate::Float64, 
    strike::Float64)
    angle = az_to_angle(strike)

    rotate_velocity(dextral_rate, extension_rate, angle)
end


function ve_vn_to_fault_slip_rate(ve::Float64, vn::Float64, strike::Float64)
    angle = az_to_angle(strike)

    rotate_velocity(ve, vn, -angle)
end


"""
    fault_oblique_merc(fault, lons, lats)
Projects a fault and other model elements from geographic coordinates into
an oblique Mercator projection defined by the fault trace.

Modified from Meade and Loveless, Blocks, `faultobliquemerc.m`.

See USGS "Map Projections - A Working Manual" p. 69 for mathematical reference.

# Arguments

# Returns

"""
function fault_oblique_merc(fault::Fault, lons::Array{Float64}, 
                            lats::Array{Float64}, R = EARTH_RAD_KM)
    
    n_stations = length(lons)

    lons = [lons; fault.trace[1,1]; fault.trace[end,1]]
    lats = [lats; fault.trace[1,2]; fault.trace[end,2]]

    lon1 = ones(n_stations + 2) .* fault.trace[1,1]
    lat1 = ones(n_stations + 2) .* fault.trace[1,2]
    lon2 = ones(n_stations + 2) .* fault.trace[end,1]
    lat2 = ones(n_stations + 2) .* fault.trace[end,2]

    # trig functions
    clat =  cosd.(lats)
    slat =  sind.(lats)
    clat1 = cosd.(lat1)
    slat1 = sind.(lat1)
    clat2 = cosd.(lat2)
    slat2 = sind.(lat2)
    clon1 = cosd.(lon1)
    slon1 = sind.(lon1)
    clon2 = cosd.(lon2)
    slon2 = sind.(lon2)

    # Pole longitude
    num = clat1 .* slat2 .* clon1 .- slat1 .* clat2 .* clon2
    den = slat1 .* clat2 .* slon2 .- clat1 .* slat2 .* slon1
    lonp = rad2deg.(atan.(num, den))

    # Pole latitutde
    latp = atand.(-cosd.(lonp .- lon1) ./ tand.(lat1))
    sp = sign.(latp)

    # choose northern hemisphere pole
    lonp[latp .< 0.] .+= 180.
    latp[latp .< 0.] .*= -1.
    # find origin longitude
    lon0 = lonp .+ 90.
    lon0[lon0 .> 180.] .-= 360.
    
    clatp = cosd.(latp)
    slatp = sind.(latp)
    dlon = lons .- lon0
    A = slatp .* slat .- clatp .* clat .* sind.(dlon)

    # Projection
    x = atan.((tand.(lats) .* clatp .+ slatp .* sind.(dlon)) ./ cosd.(dlon))
    x[latp .< 80.] = x[latp .< 80.] .- (cosd.(dlon[latp .< 80.]) .> 0.) .* pi .+ pi / 2.
    x[latp .>= 80.] = (x[latp .>= 80.] .- (cosd.(dlon[latp .>= 80.]) .< 0.) .* pi .+
                       pi / 2.)
    y = atanh.(A)

    x = -sp .* x .* R
    y = -sp .* y .* R

    (x, y)
end


"""
    build_strike_rot_matrix(strike)

Creates a rotation matrix to rotate velocities from a strike-
aligned reference frame to a cardinal direction reference frame
(and vice versa).

Called 'P_alpha' in Meade and Loveless, 2009.

# Arguments:
- `strike::Float64`: The strike of a fault.

# Returns:
- `P_strike`: 3x3 Float64 matrix.
"""
function build_strike_rot_matrix(strike::Float64)
    strike_ang = az_to_angle(strike)
    P_strike = [cos(strike_ang) -sin(strike_ang) 0.;
                sin(strike_ang)  cos(strike_ang) 0.;
                0.               0.              1.]
end

"""
    build_Pf_vert(strike)

Builds a matrix that projects differential east and north velocities
of block motions across a fault to the strike-parallel and strike-
perpendicular components.  For use with vertical faults.

Called 'P_f' in Meade and Loveless 2009.

# Arguments:
- `strike::Float64`: The strike of a fault.

# Returns:
- `P_f`: 3x3 Float64 matrix.
"""
function build_Pf_vert(strike::Float64)
    strike_ang = az_to_angle(strike)
    Pf_vert = [cos(strike_ang) -sin(strike_ang) 0.;
               0.            0.           0.;
               sin(strike_ang)  cos(strike_ang) 0.]
end


"""
    build_Pf_dip(strike, dip)

Builds a matrix that projects differential east and north velocities
of block motions across a fault to the strike-parallel and strike-
perpendicular components.  For use with non-vertical faults.

Called 'P_f' in Meade and Loveless 2009.

# Arguments:
- `strike::Float64`: The strike of a fault.
- `dip::Float64`: The dip of a fault.

# Returns:
- `P_f`: 3x3 Float64 matrix.
"""
function build_Pf_dip(strike::Float64, dip::Float64)
    strike_ang = az_to_angle(strike)
    cd = cosd(dip)

    Pf_dip = [cos(strike_ang)       -sin(strike_ang)        0.;
              sin(strike_ang) * cd   cos(strike_ang) * cd   0.;
              0.                 0.                         0.]
end


"""
    build_velocity_projection_matrix(strike, dip)

Builds a matrix that projects differential east and north velocities
of block motions across a fault to the strike-parallel and strike-
perpendicular components.

Called 'P_f' in Meade and Loveless 2009.

# Arguments:
- `strike::Float64`: The strike of a fault.
- `dip::Float64`: The dip of a fault.

# Returns:
- `P_f`: 3x3 Float64 matrix.
"""
function build_velocity_projection_matrix(strike::Float64, dip::Float64)
    # Currently we only model the horizontal components (ss, ext), so only
    # use the strike-slip formulation
    Pf = build_Pf_vert(strike)

    # if dip >= 89. # many SS faults are given 89 deg dips to have hw, fw defined
    #    Pf = build_Pf_vert(strike)
    # else
    #    Pf = build_Pf_dip(strike, dip)
    # end
    Pf
end





end # module
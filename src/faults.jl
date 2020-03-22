module Faults
export Fault, fault_to_vel, fault_slip_rate_to_ve_vn

using Parameters

using ..Oiler: VelocityVectorSphere, average_azimuth, az_to_angle, angle_difference

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
@with_kw struct Fault
    trace::Array{Float64,2}
    strike::Float64
    dip::Float64
    dip_dir::String
    extension_rate::Float64 = 0.
    extension_err::Float64 = 0.
    dextral_rate::Float64 = 0.
    dextral_err::Float64 = 0.
    lsd::Float64 = 25.
    usd::Float64 = 0.
    name::String = ""
    hw::String = ""
    fw::String = ""
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
        fix = fault.hw, mov = fault.fw, name = fault.name, ee = ee, en = en)
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



function fault_slip_rate_to_ve_vn(dextral_rate::Float64,
    extension_rate::Float64, strike::Float64)

    ang = az_to_angle(strike)

    ve = dextral_rate * cos(ang) - extension_rate * sin(ang)
    vn = dextral_rate * sin(ang) + extension_rate * cos(ang)

    (ve, vn)
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
    P_strike = [cosd(strike) -sind(strike) 0.;
                sind(strike)  cosd(strike) 0.;
                0.            0.           1.]
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
    Pf_vert = [cosd(strike) -sind(strike) 0.;
               0.            0.           0.;
               sind(strike)  cosd(strike) 0.]
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
    cd = cosd(dip)

    Pf_dip = [cosd(strike)       -sind(strike)        0.;
              sind(strike) / cd   cosd(strike) / cd   0.;
              0.                  0.                  0.]
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
    if dip == 90.
        Pf = build_Pf_vert(strike)
    else
        Pf = build_Pf_dip(strike, dip)
    end
    Pf
end





end # module
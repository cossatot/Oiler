module Faults
export Fault, fault_to_vel, fault_slip_rate_to_ve_vn, ve_vn_to_fault_slip_rate,
    fault_oblique_merc, build_strike_rot_matrix,
    build_velocity_projection_matrix, fault_to_vels

using Parameters

using Oiler
using ..Oiler: VelocityVectorSphere, average_azimuth, az_to_angle,
    angle_difference, rotate_velocity, EARTH_RAD_KM, oblique_merc, PoleCart,
    PoleSphere

const direction_map = Dict{String,Float64}("N" => 0.0,
    "NNE" => 22.5,
    "NE" => 45.0,
    "ENE" => 67.5,
    "E" => 90.0,
    "ESE" => 112.5,
    "SE" => 135.0,
    "SSE" => 157.5,
    "S" => 180.0,
    "SSW" => 202.5,
    "SW" => 225.0,
    "WSW" => 247.5,
    "W" => 270.0,
    "WNW" => 292.5,
    "NW" => 315.0,
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
    cde::Float64
    lsd::Float64
    usd::Float64
    name::String
    hw::String
    fw::String
    fid
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
- `cde::Float64 = 0.`: covariance between dextral and extension rates/errs.
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
    extension_rate::Float64=0.0,
    extension_err::Float64=0.0,
    dextral_rate::Float64=0.0,
    dextral_err::Float64=0.0,
    cde::Float64=0.0,
    lsd::Float64=10.0,
    usd::Float64=0.0,
    name::String="",
    hw::String="",
    fw::String="",
    fid="",
    check_trace=true)

    if check_trace
        trace_start = trace[1]
        if dip != 90.0
            trace = check_right_hand_rule(trace, dip_dir)
        end
        strike = average_azimuth(trace[:, 1], trace[:, 2])
        if trace[1] != trace_start
            @warn "reversing $fid trace"
        end
    else
        strike = average_azimuth(trace[:, 1], trace[:, 2])
    end

    Fault(trace, strike, dip, dip_dir,
        extension_rate, extension_err,
        dextral_rate, dextral_err, cde,
        lsd, usd, name, hw, fw, fid)
end

"""
    fault_to_vel(fault)
Converts a `Fault` object into a `VelocityVectorSphere`.  

The velocity is expressed as the horizontal velocity of the moving footwall wall
relative to the fixed hanging wall (for vertical faults, this is the moving
right-hand-side block relative to the fixed left-hand-side block when looking
from the last trace coordinate to the first trace coordinate). The geographic
location of the `VelocityVectorSphere` is the middle of the fault trace.
"""
function fault_to_vel(fault::Fault)
    ve, vn = fault_slip_rate_to_ve_vn(fault.dextral_rate, fault.extension_rate,
        fault.strike)

    ee, en, cen = fault_slip_rate_err_to_ee_en(fault.dextral_err, fault.extension_err,
        fault.strike; cde=fault.cde)

    vlon, vlat = get_midpoint(fault.trace)

    VelocityVectorSphere(lon=vlon, lat=vlat, ve=ve, vn=vn,
        fix=fault.hw, mov=fault.fw, name=fault.name, ee=ee, en=en, cen=cen,
        vel_type="fault")
end


function fault_to_vel_point(fault::Fault)
    vlon, vlat = get_midpoint(fault.trace)

    return Dict("lon" => vlon, "lat" => vlat,
        "rl" => fault.dextral_rate, "ex" => fault.extension_rate,
        "e_rl" => fault.dextral_err, "e_ex" => fault.extension_err)
end


function get_midpoint(trace::Array{Float64,2})
    fault_length = Oiler.Geom.polyline_length(trace)
    mid_pt = Oiler.Geom.sample_polyline(trace, [fault_length / 2.0])[1]
end


function get_midpoint_old(trace::Array{Float64,2})
    n_pts = size(trace)[1]

    if iseven(n_pts)
        mid = n_pts ÷ 2
        return ((trace[mid, 1] + trace[mid+1, 1]) / 2.0,
            (trace[mid, 2] + trace[mid+1, 2]) / 2.0)
    else
        mid = n_pts ÷ 2 + 1
        return (trace[mid, 1], trace[mid, 2])
    end
end


function check_right_hand_rule(trace::Array{Float64,2}, dip_dir::String;
    reverse_angle_threshold::Float64=90.0)
    # Modified from the OQ-MBTK tools, (c) Global Earthquake Model Foundation

    strike = average_azimuth(trace[:, 1], trace[:, 2])

    trace_dip_trend = strike + 90.0

    fault_dip_trend = direction_map[dip_dir]

    trend_angle_difference = angle_difference(trace_dip_trend, fault_dip_trend)

    # TODO: warn if trend_angle_difference < 15 (or some other low threshold)
    if trend_angle_difference > reverse_angle_threshold
        trace = reverse(trace, dims=1)
    end
    trace
end


function fault_slip_rate_to_ve_vn(dextral_rate::Float64, extension_rate::Float64,
    strike::Float64)
    angle = az_to_angle(strike)

    rotate_velocity(dextral_rate, extension_rate, angle)
end

function fault_slip_rate_err_to_ee_en(dextral_err::Float64, extension_err::Float64,
    strike::Float64; cde::Float64=0.0)
    angle = az_to_angle(strike)

    Oiler.Geom.rotate_velocity_err(dextral_err, extension_err, angle; cov=cde)
end


function ve_vn_to_fault_slip_rate(ve::Float64, vn::Float64, strike::Float64)
    angle = az_to_angle(strike)

    rotate_velocity(ve, vn, -angle)
end


function ee_en_to_fault_slip_rate_err(ee::Float64, en::Float64, strike::Float64;
    cen::Float64=0.0)
    angle = az_to_angle(strike)

    Oiler.Geom.rotate_velocity_err(ee, en, -angle; cov=cen)
end


function get_fault_slip_rate_from_pole(fault::Oiler.Faults.Fault, pole::Oiler.PoleCart;
    lon=nothing, lat=nothing)

    if isnothing(lon)
        fault_mid = get_midpoint(fault.trace)
        lon, lat = fault_mid[1], fault_mid[2]
    end

    # fault velocities should be relative to hanging wall (hw fixed)
    if (fault.fw == pole.mov) & (fault.hw == pole.fix)
        pole_use = pole
    elseif (fault.fw == pole.fix) & (fault.hw == pole.mov)
        pole_use = -pole
    else
        fw = fault.fw
        hw = fault.hw
        mov = pole.mov
        fix = pole.fix
        warn_msg = "fault ($hw, $fw) and pole ($fix, $mov) do not match, but continuing..."
        @warn warn_msg
        pole_use = pole
    end

    pred_vel = Oiler.BlockRotations.predict_block_vel(lon, lat, pole_use)
    v_rl, v_ex = ve_vn_to_fault_slip_rate(pred_vel.ve, pred_vel.vn, fault.strike)

    if (pred_vel.ee == 0.0) & (pred_vel.en == 0.0)
        e_rl, e_ex, cde = 0.0, 0.0, 0.0
    else
        e_rl, e_ex, cde = ee_en_to_fault_slip_rate_err(pred_vel.ee, pred_vel.en,
            fault.strike; cen=pred_vel.cen)
    end
    v_rl, v_ex, e_rl, e_ex, cde
end


"""
    fault_to_vels

Makes multiple velocities from a fault. This is to weight the velocity inversion
by the length of faults, and to provide some damping against excessive
local rotations induced by too few local constraints.
"""
function fault_to_vels(fault::Fault; simp_dist::Float64=1.0)
    simp_trace = Oiler.Geom.simplify_polyline(fault.trace, simp_dist)
    fault_length = Oiler.Geom.polyline_length(simp_trace)

    if fault_length < 10.0
        n_segs = 1
    elseif (10.0 <= fault_length) & (fault_length < 30.0)
        n_segs = 2
    elseif (30.0 <= fault_length) & (fault_length < 60.0)
        n_segs = 3
    else
        #n_segs = Int(floor(150. / 50.))
        n_segs = Int(floor(fault_length / 20.0))
    end

    new_traces = Oiler.Geom.break_polyline_equal(simp_trace, n_segs)

    vels = [fault_to_vel(
        Fault(; trace=tr, dip=fault.dip, dip_dir=fault.dip_dir,
            extension_rate=fault.extension_rate,
            extension_err=fault.extension_err,
            dextral_rate=fault.dextral_rate, dextral_err=fault.dextral_err,
            lsd=fault.lsd, usd=fault.usd, name=fault.name, hw=fault.hw,
            fw=fault.fw))
            for tr in new_traces]
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
    lats::Array{Float64})

    n_stations = length(lons)

    lons = [lons; fault.trace[1, 1]; fault.trace[end, 1]]
    lats = [lats; fault.trace[1, 2]; fault.trace[end, 2]]

    # lon1 = ones(n_stations + 2) .* fault.trace[1,1]
    # lat1 = ones(n_stations + 2) .* fault.trace[1,2]
    # lon2 = ones(n_stations + 2) .* fault.trace[end,1]
    # lat2 = ones(n_stations + 2) .* fault.trace[end,2]

    lon1 = fault.trace[1, 1]
    lat1 = fault.trace[1, 2]
    lon2 = fault.trace[end, 1]
    lat2 = fault.trace[end, 2]

    oblique_merc(lons, lats, lon1, lat1, lon2, lat2)

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
    P_strike = [cos(strike_ang) -sin(strike_ang) 0.0
        sin(strike_ang) cos(strike_ang) 0.0
        0.0 0.0 1.0]
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
    Pf_vert = [cos(-strike_ang) -sin(-strike_ang) 0.0
        0.0 0.0 0.0
        0.0 0.0 0.0]
    # sin(-strike_ang)  cos(-strike_ang) 0.]
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
    # cd = cosd(dip)
    cd = 1.0  # reasoning: longer-term, all convergence/extension goes down-dip

    Pf_dip = [cos(-strike_ang) -sin(-strike_ang) 0.0
        sin(-strike_ang)/cd cos(-strike_ang)/cd 0.0
        0.0 0.0 0.0]
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
    # Pf = build_Pf_vert(strike)

    if dip >= 89.0 # many SS faults are given 89 deg dips to have hw, fw defined
        Pf = build_Pf_vert(strike)
    else
        Pf = build_Pf_dip(strike, dip)
    end
    Pf
end


function project_fault_trace(fault)
    hdist = fault.lsd / tand(fault.dip)
    az = fault.strike + 90.0

    new_trace = zeros(size(fault.trace))

    for i = 1:size(fault.trace, 1)
        new_trace[i, :] .= Oiler.Geom.terminal_coords_from_bearing_dist(
            fault.trace[i, 1], fault.trace[i, 2], az, hdist)
    end
    new_trace
end


end # module

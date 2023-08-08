module Geom

export azimuth, gc_distance, average_azimuth, az_to_angle, angle_to_az,
    angle_difference, rotate_velocity, rotate_xy_vec, oblique_merc
import Statistics: mean
#import Proj: CRS, Transformation

using Oiler
using ..Oiler: EARTH_RAD_KM
using LinearAlgebra

using Proj: CRS, Transformation

struct LineString
    coords::Array{Float64,2}
end

struct Polygon
    coords::Array{Float64,2}
end

struct Point
    coords::Array{Float64,2}
end


function azimuth(lon1::Float64, lat1::Float64,
    lon2::Float64, lat2::Float64)

    dlon = lon2 - lon1
    y = sind(dlon) * cosd(lat2)
    x = cosd(lat1) * sind(lat2) - sind(lat1) * cosd(lat2) * cosd(dlon)

    azimuth = atand(y, x)
end


function gc_distance(lon1::Float64, lat1::Float64,
    lon2::Float64, lat2::Float64, R=EARTH_RAD_KM)

    lon1 = deg2rad(lon1)
    lat1 = deg2rad(lat1)
    lon2 = deg2rad(lon2)
    lat2 = deg2rad(lat2)

    distance = asin(sqrt(sin((lat1 - lat2) / 2.0)^2.0
                         +
                         cos(lat1) * cos(lat2)
                         * sin((lon1 - lon2) / 2.0)^2.0))

    distance * 2.0 * R
end


function angle_fixer_deg(ang)
    while ang < 0
        ang += 360.0
    end

    while ang > 360.0
        ang -= 360.0
    end
    ang
end


function average_azimuth(lons::Array{Float64}, lats::Array{Float64})

    if length(lons) == 2
        az = azimuth(lons[1], lats[1], lons[2], lats[2])
    else

        azs = [azimuth(lons[i-1], lats[i-1], lons[i], lats[i]) for i =
            2:length(lats)]

        azs = [deg2rad(az) for az in azs]

        dists = [gc_distance(lons[i-1], lats[i-1], lons[i], lats[i]) for i =
            2:length(lats)]

        avg_x = mean(dists .* [sin(az) for az in azs])
        avg_y = mean(dists .* [cos(az) for az in azs])

        az = rad2deg(atan(avg_x, avg_y))
    end

    az = angle_fixer_deg(az)

    az
end


function average_azimuth(linestring::LineString)
    
    coords = linestring.coords
    
    return Oiler.Geom.average_azimuth(coords[:, 1], coords[:, 2])
end


function check_duplicates(arr, name="")
    for i = 1:size(arr, 1)-1
        if arr[i] == arr[i+1]
            @warn "duplicated consecutive values in $name"
        end
    end
end


function az_to_angle(az::Float64)
    deg2rad(90.0 - az)
end

function angle_to_az(angle::Float64)
    az = -(rad2deg(angle) - 90.0)
    az = angle_fixer_deg(az)
    az
end


"""
    angle_difference(trend_1, trend_2)
Calculates the difference between to angles, azimuths or trends, in degrees.
"""
function angle_difference(trend_1::Float64, trend_2::Float64; return_abs::Bool=true,
    unit="degrees")

    if unit == "radians"
        wrap_min = -pi
        wrap_max = 2 * pi
    elseif unit == "degrees"
        wrap_min = -180.0
        wrap_max = 360.0
    end

    difference = trend_2 - trend_1

    while difference < wrap_min
        difference += wrap_max
    end
    while difference > -wrap_min
        difference -= wrap_max
    end

    if return_abs == true
        difference = abs(difference)
    end
    difference
end


function colatitude(latitude)
    return angle_difference(90., latitude)
end


function angular_mean_degrees(angles)
    mean_angle = rad2deg(atan(mean(sind.(angles)), mean(cosd.(angles))))
end

function angular_mean_radians(angles)
    mean_angle = atan(mean(sin.(angles)), mean(cos.(angles)))
end

function angular_mean(angles; unit="radians")
    if unit == "radians"
        mean_angle = angular_mean_radians(angles)
    elseif unit == "degrees"
        mean_angle = angular_mean_degrees(angles)
    else
        throw(ArgumentError("Unit needs to be radians or degrees"))
    end
    mean_angle
end


function angular_std(angles; unit="radians")
    mean = angular_mean(angles; unit=unit)

    diffs = [angle_difference(angle, mean; return_abs=false, unit=unit)
             for angle in angles]

    var = sum(diffs .^ 2) / (length(diffs) - 1)
    std = sqrt(var)
end


function point_sphere_to_cart(lon, lat, depth; R=EARTH_RAD_KM)
    r = R + depth

    x = r * cosd(lat) * cosd(lon)
    y = r * cosd(lat) * sind(lon)
    z = r * sind(lat)

    [x, y, z]
end


function point_sphere_to_cart(pos; R=EARTH_RAD_KM)
    lon = pos[1]
    lat = pos[2]
    depth = pos[3]

    r = R + depth

    x = r * cosd(lat) * cosd(lon)
    y = r * cosd(lat) * sind(lon)
    z = r * sind(lat)

    [x, y, z]
end


function point_cart_to_sphere(pos; R=EARTH_RAD_KM)
    x = pos[1]
    y = pos[2]
    z = pos[3]

    point_cart_to_sphere(x, y, z; R=R)
end


function point_cart_to_sphere(x, y, z; R=EARTH_RAD_KM)
    r = sqrt(x^2 + y^2 + z^2)
    depth = r - R

    lon = rad2deg(atan(y, x))
    lat = rad2deg(asin(z / r))

    [lon, lat, depth]
end


function rotate_velocity(vx::Float64, vy::Float64, angle::Float64)
    Vp = [cos(angle) -sin(angle); sin(angle) cos(angle)] * [vx; vy]
end


function rotate_velocity_err(ex, ey, angle; cov=0.0)

    R = [cos(angle) -sin(angle); sin(angle) cos(angle)]

    M = [ex^2 cov; cov ey^2]

    var = R * M * R'
    err = [sqrt(var[1, 1]), sqrt(var[2, 2]), var[1, 2]]
end


function rotate_velocity_w_err(vx::Float64, vy::Float64, angle::Float64,
    ex=0.0, ey=0.0; cov=0.0)
    Vp = [cos(angle) -sin(angle); sin(angle) cos(angle)] * [vx; vy]
    Ep = rotate_velocity_err(ex, ey, angle; cov=cov)
    (Vp, Ep)
end


"""
    rotate_xy_vec

Rotates x and y vectors based on a rotation angle in radians.
"""
function rotate_xy_vec(x, y, alpha_rot)
    xp = cos(alpha_rot) .* x - sin(alpha_rot) .* y
    yp = sin(alpha_rot) .* x + cos(alpha_rot) .* y
    return xp, yp
end


function rotate_xy_vec_to_magnitude(x, y; x_err=0.0, y_err=0.0, cov=0.0)
    angle = atan(-y, x)
    (vec_rot, err_rot) = rotate_velocity_w_err(x, y, angle, x_err, y_err; cov=cov)
    (vec_rot[1], err_rot[1])
end

function rotate_xy_vec_to_magnitude(x::Array{Float64}, y::Array{Float64}; 
        x_err=nothing, y_err=nothing, cov=nothing)
    angle = atan.(-y, x)

    if isnothing(x_err)
        x_err = zeros(size(x))
    end
    if isnothing(y_err)
        y_err = zeros(size(x))
    end
    if isnothing(cov)
        cov = zeros(size(x))
    end

    mags = zeros(size(x))
    errs = zeros(size(x))

    for (i, x_) in enumerate(x)

        (vec_rot, err_rot) = rotate_velocity_w_err(x_, y[i], angle[i], 
                                                   x_err[i], y_err[i]; 
                                                   cov=cov[i])
        mags[i] = vec_rot[1]
        errs[i] = err_rot[1]

    end
    mags, errs
end


function get_oblique_merc(lon1, lat1, lon2, lat2)
    # correction for perfectly horizontal lines or lat1 at zero
    if (abs(lat1 - lat2) < 2e-3) || (lat1 == 0.0)
        lat1 = lat1 + 2e-3
    end

    init_str = "+proj=omerc +lat_1=$lat1 +lon_1=$lon1 +lat_2=$lat2 +lon_2=$lon2 +ellps=WGS84"
end


function oblique_merc(lons, lats, lon1, lat1, lon2, lat2)
    wgs84 = "+proj=longlat +datum=WGS84 +nodefs"
    omerc = get_oblique_merc(lon1, lat1, lon2, lat2)
    trans = Transformation(wgs84, omerc; always_xy=true)


    x = zeros(size(lons))
    y = zeros(size(lons))

    for i in eachindex(lons)
        @inbounds x[i], y[i] = trans(lons[i], lats[i])
    end

    (x, y)
end


function bearing(final_lon, final_lat, start_lon, start_lat)

    y = sind(final_lon - start_lon) * cosd(final_lat)
    x = cosd(start_lat) * sind(final_lat) - sind(start_lat) * cosd(final_lat) * cosd(final_lon-start_lon)

    atand(y,x)
end


function azimuthal_equidistant_proj(lons, lats, clon, clat)


    dists = zeros(size(lons))
    angles = zeros(size(lons))

    #dists = [gc_distance(lon, lats[i], clon, clat) for (i, lon) in enumerate(lons)]
    for (i, lon) in enumerate(lons)
        @inbounds dists[i] += gc_distance(lon, lats[i], clon, clat) * 1000.
        @inbounds angles[i] += az_to_angle(bearing(lon, lats[i], clon, clat))
    end

    xs = dists .* cos.(angles)
    ys = dists .* sin.(angles)

    xs, ys
end


function azimuthal_equidistant_proj(lons, lats, lon1, lat1, lon2, lat2)
    clon, clat = sample_polyline([lon1 lat1; lon2 lat2], 
                                 [gc_distance(lon1, lat1, lon2, lat2)/2.])[1]

    azimuthal_equidistant_proj(lons, lats, clon, clat)
end


function calculate_central_meridian(point1, point2, point3)
    longitudes_rad = deg2rad.([point1[1], point2[1], point3[1]])

    mean_sin = mean(sin.(longitudes_rad))
    mean_cos = mean(cos.(longitudes_rad))
    
    central_meridian = rad2deg(atan(mean_sin, mean_cos))
    
    # Normalize to the range [-180, 180)
    central_meridian = (central_meridian + 180) % 360 - 180

    return central_meridian
end


function transverse_mercator_projection(lon, lat, central_meridian)
    R = EARTH_RAD_KM * 1000.
    phi = deg2rad(lat)
    lambda_ = deg2rad(lon)
    lambda_0 = deg2rad(central_meridian)

    x = R * (lambda_ - lambda_0) * cos(phi)
    y = R * log(tan(π/4 + phi/2))

    return x, y
end


"""
    terminal_coords_from_bearing_dist

Finds the coordinates for a point that is some distance
away from another point at a specified (initial) bearing.

Formula from https://www.movable-type.co.uk/scripts/latlong.html
"""
function terminal_coords_from_bearing_dist(lon1, lat1, bearing, dist)
    ang_dist = dist / EARTH_RAD_KM

    lat2 = asind(sind(lat1) * cos(ang_dist) + cosd(lat1) * sin(ang_dist) * cosd(bearing))
    lon2 = lon1 + atand(sind(bearing) * sin(ang_dist) * cosd(lat1),
        cos(ang_dist) - sind(lat1) * sind(lat2))
    lon2, lat2
end


function polyline_seg_lengths(polyline::Array{Float64,2})
    n_segs = size(polyline, 1) - 1
    seg_lengths = Array{Float64}(undef, n_segs)

    for i = 1:n_segs
        seg_lengths[i] = gc_distance(polyline[i, 1], polyline[i, 2],
            polyline[i+1, 1], polyline[i+1, 2])
    end
    seg_lengths
end


"""
    polyline_length(polyline)

Calculates the length of a polyline of (longitude, latitude) coordinates
as the sum of great-circle distances between the points.  Length is returned
in kilometers.
"""
function polyline_length(polyline::Array{Float64,2})
    seg_lengths = polyline_seg_lengths(polyline)
    sum(seg_lengths)
end


"""
    sample_polyline(polyline, dists)

Returns the coordinates of points along a polyline at specified distances
from the start point. The line should be in longitude, latitude (WGS84) format
and the distances should be in km.  All distance calculations are the great
circle distance.

# Arguments
- `polyline::Array{Float64,2}``: An array of coordinates in (lon, lat) form.
- `dists`: Some collection or iterable of distances from the start point. These
    should be positive floats.

# Returns
Coordinates for each positive distance that is less than the total length
of the polyline. Note that no coordinates will be returned for distances
that are negative or greater than the total length of the line, so the
number of returned coordinates may be less than the number of distances passed
as the argument.
"""
function sample_polyline(polyline::Array{Float64,2}, dists)
    seg_lengths = polyline_seg_lengths(polyline)
    n_segs = length(seg_lengths)
    cum_lengths = vcat(0.0, cumsum(seg_lengths))

    new_pts = Array{Float64,2}[]

    for dist in dists[dists.>0]
        last_smaller_idx = size(cum_lengths[cum_lengths.<dist], 1)
        if last_smaller_idx <= n_segs
            start_lon, start_lat = polyline[last_smaller_idx, :]
            end_lon, end_lat = polyline[last_smaller_idx+1, :]
            seg_az = azimuth(start_lon, start_lat, end_lon, end_lat)
            remain_dist = dist - cum_lengths[last_smaller_idx]

            new_lon, new_lat = terminal_coords_from_bearing_dist(start_lon,
                start_lat, seg_az, remain_dist)

            push!(new_pts, [new_lon new_lat])
            # elseif ??
        end
    end
    vcat(new_pts)
end


"""
Simplifies a polyline based on the Ramer-Douglas-Peucker algorithm.

Code based on Fabian Hirschmann's RDP Python code.
"""
function simplify_polyline(polyline, min_dist)
    dmax = 0.0
    n = size(polyline)[1]
    index = n * 1

    for i = 1:n
        d = point_line_distance(polyline[i, :], polyline[1, :], polyline[n, :])

        if d > dmax
            index = i
            dmax = d
        end
    end

    if dmax > min_dist
        r1 = simplify_polyline(polyline[1:index, :], min_dist)
        r2 = simplify_polyline(polyline[index:n, :], min_dist)

        simp_line = [r1[1:end-1, :]; r2]
    else
        simp_line = [polyline[1:1, :]; polyline[n:n, :]]
    end
    simp_line
end


function two_cross(vec1, vec2)
    vec1[1] * vec2[2] - vec1[2] * vec2[1]
end

function point_line_distance(point, line_start, line_end)
    if (line_start[1] == line_end[1])
        if (line_start[2] == line_end[2])
            return norm(point .- line_start)
        end
    end

    top = abs(norm(two_cross(line_end .- line_start, line_start .- point)))
    bottom = norm(line_end .- line_start)
    return top / bottom
end



function break_polyline_equal(polyline, n_segs)
    segment_lengths = polyline_seg_lengths(polyline)
    vertex_cum_dists = vcat(0.0, cumsum(segment_lengths))
    poly_length = vertex_cum_dists[end]

    new_polylines = []
    if n_segs > 1
        seg_length = poly_length / n_segs
        cum_lengths = seg_length .* collect(1:n_segs)
        break_dists = cum_lengths[1:end-1]

        break_coords = sample_polyline(polyline, break_dists)
        break_idxs = [length(vertex_cum_dists[vertex_cum_dists.<break_dist])
                      for break_dist in break_dists]

        for n_seg = 1:n_segs
            if n_seg == 1
                start_idx = 1
                stop_idx = break_idxs[1]
                seg_trace = vcat(polyline[start_idx:stop_idx, :],
                    break_coords[1])
            elseif (n_seg > 1) & (n_seg < n_segs)
                start_idx = break_idxs[n_seg-1] + 1
                stop_idx = break_idxs[n_seg]
                seg_trace = vcat(break_coords[n_seg-1],
                    polyline[start_idx:stop_idx, :],
                    break_coords[n_seg])
            else
                start_idx = break_idxs[n_seg-1] + 1
                seg_trace = vcat(break_coords[n_seg-1],
                    polyline[start_idx:end, :])
            end # if
            push!(new_polylines, seg_trace)
        end # for n_seg
    else
        push!(new_polylines, polyline)
    end
    new_polylines
end

"""
   returns 1 for CCW, and -1 for CW 
"""
function check_winding_order(coords::Array{Float64,2})

    function fun(p1, p2)
        (p2[1] - p1[1]) * (p2[2] + p1[2])
    end

    Int(sign(sum([fun(coords[i, :], coords[i-1, :]) for i = 2:size(coords, 1)])))
end


"""
   returns 1 for CCW, and -1 for CW 
"""
function check_winding_order(coords)

    function fun(p1, p2)
        (p2[1] - p1[1]) * (p2[2] + p1[2])
    end

    Int(sign(sum([fun(coords[i], coords[i-1]) for i = 2:size(coords, 1)])))
end


function strike_dip_from_3_pts(pt1, pt2, pt3)
    A = [pt1[1] pt1[2] 1.0
        pt2[1] pt2[2] 1.0
        pt3[1] pt3[2] 1.0]

    zvec = [-pt1[3], -pt2[3], -pt3[3]]

    mx, my, z0 = A \ zvec

    dip = atand(sqrt(mx^2 + my^2))
    dip_dir = atan(my, mx)
    strike_angle = dip_dir + pi / 2.0
    # strike_angle = dip_dir
    strike = angle_to_az(strike_angle)

    strike, dip

end


function get_coord_vecs(coords::Array{Point})
    lons = collect(p.coords[1] for p in coords)
    lats = collect(p.coords[2] for p in coords)

    lons, lats
end


function polygon_center(pts)
    centroid = [0.0 0.0]

    n = size(pts, 1)
    signed_area = 0.0

    for i in 1:n
        ii = i + 1
        if ii > n
            ii = 1
        end
        x0, y0 = pts[i, :]
        x1, y1 = pts[ii, :]

        A = (x0 * y1) - (x1 * y0)
        signed_area += A

        centroid[1] += (x0 + x1) * A
        centroid[2] += (y0 + y1) * A
    end

    centroid ./= 3 * signed_area
end


function get_polygon_centroid(poly::Polygon; epsg=102016, algorithm=polygon_center)

    poly_pts = poly.coords

    if epsg != 4326
        trans = Oiler.IO.make_trans_from_wgs84(epsg)
        poly_zeros = zeros(size(poly_pts))
        for (i, row) in enumerate(eachrow(poly_pts))
            poly_zeros[i, :] = collect(trans(row))
        end
        poly_pts = poly_zeros
    end

    centroid = algorithm(poly_pts)

    if epsg != 4326
        back_trans = Oiler.IO.make_trans_to_wgs84(epsg)
        bc = back_trans(centroid)
        centroid = [bc[1] bc[2]]
    end

    Point(centroid)
end

end # module

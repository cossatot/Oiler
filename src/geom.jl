module Geom

export azimuth, gc_distance, average_azimuth, az_to_angle, angle_to_az,
    angle_difference, rotate_velocity, rotate_xy_vec, oblique_merc

import Statistics:mean
import Proj4: Projection, transform

using ..Oiler:EARTH_RAD_KM
using LinearAlgebra

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

    distance = asin(sqrt(sin((lat1 - lat2) / 2.)^2.
        + cos(lat1) * cos(lat2)
        * sin((lon1 - lon2) / 2.)^2.))

    distance * 2. * R 
end


function average_azimuth(lons::Array{Float64}, lats::Array{Float64})

    if length(lons) == 2
        az = azimuth(lons[1], lats[1], lons[2], lats[2])
    else
        
        azs = [azimuth(lons[1], lats[1], lons[2], lats[2]) for i in
            2:length(lats)]
        
        azs = [deg2rad(az) for az in azs]

        dists = [gc_distance(lons[1], lats[1], lons[i], lats[i]) for i in
            2:length(lats)]

        avg_x = mean(dists .* [sin(az) for az in azs])
        avg_y = mean(dists .* [cos(az) for az in azs])

        az = rad2deg(atan(avg_x, avg_y))
        if az < 0.
            az += 360.
        end
    end
    az
end
        

function az_to_angle(az::Float64)
    deg2rad(90. - az)
end

function angle_to_az(angle::Float64)
    -(rad2deg(angle) - 90.)
end

"""
    angle_difference(trend_1, trend_2)
Calculates the difference between to angles, azimuths or trends, in degrees.
"""
function angle_difference(trend_1::Float64, trend_2::Float64; return_abs::Bool=true)

    difference = trend_2 - trend_1

    while difference < -180.
        difference += 360.
    end
    while difference > 180.
        difference -= 360.
    end

    if return_abs == true
        difference = abs(difference)
    end
    difference
end


function rotate_velocity(vx::Float64, vy::Float64, angle::Float64)
    Vp = [cos(angle) -sin(angle); sin(angle) cos(angle)] * [vx; vy]
end


function rotate_velocity_err(ex, ey, angle; cov=0.)
    
    R = [cos(angle) -sin(angle); sin(angle) cos(angle)]

    M = [ex^2 cov; cov ey^2]
    
    var = R * M * R'
    err = [sqrt(var[1,1]), sqrt(var[2,2]), var[1,2]]
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


function oblique_merc(lons, lats, lon1, lat1, lon2, lat2)
    # correction for perfectly horizontal lines
    if lat1 == lat2
        lat1 = lat1 + 1e-4
    end

    wgs84 = Projection("+proj=longlat +datum=WGS84 +nodefs")
    init_str = "+proj=omerc +lat_1=$lat1 +lon_1=$lon1 +lat_2=$lat2 +lon_2=$lon2"
    omerc = Projection(init_str)

    xy = [transform(wgs84, omerc, [lon, lats[i]]) for (i, lon) in enumerate(lons)]
    x = [c[1] for c in xy]
    y = [c[2] for c in xy]
    
    (x, y)
end


"""
    terminal_coords_from_bearing_dist

Finds the coordinates for a point that is some distance
away from another point at a specified (initial) bearing.

Formula from https://www.movable-type.co.uk/scripts/latlong.html
"""
function terminal_coords_from_bearing_dist(lon1, lat1, bearing, dist)
    ang_dist = dist /  EARTH_RAD_KM

    lat2 = asind(sind(lat1) * cos(ang_dist) + cosd(lat1) * sin(ang_dist) * cosd(bearing))
    lon2 = lon1 + atand( sind(bearing) * sin(ang_dist) * cosd(lat1),
                         cos(ang_dist) - sind(lat1) * sind(lat2))
    lon2, lat2
end


function polyline_seg_lengths(polyline::Array{Float64,2})
    n_segs = size(polyline, 1) - 1
    seg_lengths = Array{Float64}(undef, n_segs)
    
    for i in 1:n_segs
        seg_lengths[i] = gc_distance(polyline[i,1], polyline[i,2],
            polyline[i + 1,1], polyline[i + 1,2])
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
    cum_lengths = vcat(0., cumsum(seg_lengths))

    new_pts = Array{Float64,2}[]
    
    for dist in dists[dists .> 0]
        last_smaller_idx = size(cum_lengths[cum_lengths .< dist], 1) 
        if last_smaller_idx <= n_segs
            start_lon, start_lat = polyline[last_smaller_idx, :]
            end_lon, end_lat = polyline[last_smaller_idx + 1, :]
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

    for i in 1:n
        d = point_line_distance(polyline[i,:], polyline[1,:], polyline[n,:])

        if d > dmax
            index = i
            dmax = d
        end
    end

    if dmax > min_dist
        r1 = simplify_polyline(polyline[1:index,:], min_dist)
        r2 = simplify_polyline(polyline[index:n,:], min_dist)

        simp_line = [r1[1:end - 1,:]; r2]
    else
        simp_line = [polyline[1:1,:]; polyline[n:n,:]]
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
    vertex_cum_dists = vcat(0., cumsum(segment_lengths))
    poly_length = vertex_cum_dists[end]

    new_polylines = []
    if n_segs > 1
        seg_length = poly_length / n_segs
        cum_lengths = seg_length .* collect(1:n_segs)
        break_dists = cum_lengths[1:end - 1]

        break_coords = sample_polyline(polyline, break_dists)
        break_idxs = [length(vertex_cum_dists[vertex_cum_dists .< break_dist])
            for break_dist in break_dists]

        for n_seg in 1:n_segs
            if n_seg == 1
                start_idx = 1
                stop_idx = break_idxs[1]
                seg_trace = vcat(polyline[start_idx:stop_idx, :], 
                                 break_coords[1])
            elseif (n_seg > 1) & (n_seg < n_segs)
                start_idx = break_idxs[n_seg - 1] + 1
                stop_idx = break_idxs[n_seg]
                seg_trace = vcat(break_coords[n_seg - 1], 
                                 polyline[start_idx:stop_idx, :],
                                 break_coords[n_seg])
            else
                start_idx = break_idxs[n_seg - 1] + 1
                seg_trace = vcat(break_coords[n_seg - 1],
                                 polyline[start_idx:end,:])
            end # if
            push!(new_polylines, seg_trace)
        end # for n_seg
    else
        push!(new_polylines, polyline)
    end
    new_polylines
end


function check_winding_order(coords::Array{Float64,2})

    function fun(p1, p2)
        (p2[1] - p1[1]) * (p2[2] + p1[2])
    end

    Int(sign(sum([fun(coords[i,:], coords[i - 1,:]) for i in 2:size(coords, 1)])))
end


function check_winding_order(coords)

    function fun(p1, p2)
        (p2[1] - p1[1]) * (p2[2] + p1[2])
    end

    Int(sign(sum([fun(coords[i], coords[i - 1]) for i in 2:size(coords, 1)])))
end



function point_in_spherical_poly(polyline)


end





end # module

module Geom

export azimuth, gc_distance, average_azimuth, az_to_angle, angle_to_az,
    angle_difference, rotate_velocity, rotate_xy_vec, oblique_merc

import Statistics: mean
import Proj4: Projection, transform

using ..Oiler: EARTH_RAD_KM
using LinearAlgebra

function azimuth(lond1::Float64, latd1::Float64,
                 lond2::Float64, latd2::Float64)

    dlon = lond2 - lond1
    y = sind(dlon) * cosd(latd2)
    x = cosd(latd1) * sind(latd2) - sind(latd1) * cosd(latd2) * cosd(dlon)

    azimuth = atand(y, x)
end


function gc_distance(lond1::Float64, latd1::Float64,
                 lond2::Float64, latd2::Float64, R = EARTH_RAD_KM)

    lon1 = deg2rad(lond1)
    lat1 = deg2rad(latd1)
    lon2 = deg2rad(lond2)
    lat2 = deg2rad(latd2)

    distance = asin(sqrt(sin((lat1 - lat2) / 2.)^2.
        + cos(lat1) * cos(lat2)
        * sin((lon1 - lon2) / 2.)^2.))

    distance * 2. * R 
end


function average_azimuth(londs::Array{Float64}, latds::Array{Float64})

    if length(londs) == 2
        az = azimuth(londs[1], latds[1], londs[2], latds[2])
    else
        
        azs = [azimuth(londs[1], latds[1], londs[2], latds[2]) for i in
            2:length(latds)]
        
        azs = [deg2rad(az) for az in azs]

        dists = [gc_distance(londs[1], latds[1], londs[i], latds[i]) for i in
            2:length(latds)]

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
function angle_difference(trend_1::Float64, trend_2::Float64; return_abs::Bool = true)

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

    (Vp[1], Vp[2])
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


function xyz2enumat(G, az)

# az = 

end


function oblique_merc(lons, lats, lon1, lat1, lon2, lat2)
    wgs84 = Projection("+proj=longlat +datum=WGS84 +nodefs")
    init_str = "+proj=omerc +lat_1=$lat1 +lon_1=$lon1 +lat_2=$lat2 +lon_2=$lon2"
    omerc = Projection(init_str)

    xy = [transform(wgs84, omerc, [lon, lats[i]]) for (i, lon) in enumerate(lons)]
    x = [c[1] for c in xy]
    y = [c[2] for c in xy]

    #x .*= 1000.
    #y .*= 1000.
    (x,y)
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

        simp_line = [r1[1:end-1,:]; r2]
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


end # module

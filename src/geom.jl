module Geom

export azimuth, gc_distance, average_azimuth, az_to_angle, angle_to_az,
    angle_difference, rotate_velocity, rotate_xy_vec, oblique_merc

import Statistics: mean
import Proj4: Projection, transform

using ..Oiler: EARTH_RAD_KM

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

    #x, y = transform(wgs84, omerc, lons, lats)

    xy = [transform(wgs84, omerc, [lon, lats[i]]) for (i, lon) in enumerate(lons)]
    x = [c[1] for c in xy]
    y = [c[2] for c in xy]

    #x .*= 1000.
    #y .*= 1000.
    (x,y)
end


function oblique_merc_bad(lons, lats, lon1, lat1, lon2, lat2, R = EARTH_RAD_KM)

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
    lonp = atand.(num, den)

    # Pole latitutde
    latp = atand.(-cosd.(lonp .- lon1) ./ tand.(lat1))
    sp = sign.(latp)

    # choose northern hemisphere pole
    #lonp[latp .< 0.] .+= 180.
    #latp[latp .< 0.] .*= -1.

    if latp < 0.
        lonp += 180.
        latp *= -1.
    end

    # find origin longitude
    lon0 = lonp .+ 90.
    #lon0[lon0 .> 180.] .-= 360.
    if lon0 > 180.
        lon0 -= 360.
    end
    
    clatp = cosd.(latp)
    slatp = sind.(latp)
    dlon = lons .- lon0
    A = slatp .* slat .- clatp .* clat .* sind.(dlon)

    # Projection
    x = atan.((tand.(lats) .* clatp .+ slatp .* sind.(dlon)) ./ cosd.(dlon))
    #x[latp .< 80.] = x[latp .< 80.] .- (cosd.(dlon[latp .< 80.]) .> 0.) .* pi .+ pi / 2.
    #x[latp .>= 80.] = (x[latp .>= 80.] .- (cosd.(dlon[latp .>= 80.]) .< 0.) .* pi .+
    #                   pi / 2.)

    #if latp < 80
    #    x .- (cosd.(dlon) .> 0.) .* pi .+ pi / 2.
    #else
    #    x = (x .- (cosd.(dlon) .< 0.) .* pi .+ pi / 2.)
    #end

    y = rad2deg.(atanh.(A))
    # Y DOES NOT AGREE WITH PROJ4

    x = -sp .* x .* R
    y = -sp .* y .* R

    (x, y)
end


end # module

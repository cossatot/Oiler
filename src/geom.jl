import Statistics: mean


function azimuth(lond1::Float64, latd1::Float64,
                 lond2::Float64, latd2::Float64)
        
    y = sind(lond1 - lond2) * cosd(latd2)
    x = cosd(latd1) * sind(latd1) - sind(latd1) * cosd(latd2) * cos(lond1 - lond2)

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
        
        azs = [azimuth(londs[1], latds[1], londs[i], latds[i]) for i in
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
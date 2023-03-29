module Strain

using ..Oiler
import Oiler.Constants.EARTH_RAD_MM


function make_strain_rate_equations(e_vel, n_vel, lon, lat, c_lon, c_lat)
    
    rlon = deg2rad(lon)
    rlat = deg2rad(Oiler.Geom.colatitude(lat)) # colatitude
    rc_lon = deg2rad(c_lon)
    rc_lat = deg2rad(Oiler.Geom.colatitude(c_lat)) # colatitude
    
    lhs = zeros(2,3)
    rhs = [e_vel; n_vel]
    
    lhs[1,1] = EARTH_RAD_MM * (rlon - rc_lon) * sin(rc_lat)  #E_lon_lon
    lhs[1,2] = EARTH_RAD_MM * (rlat - rc_lat)                #E_lon_lat
    
    lhs[2,2] = EARTH_RAD_MM * (rlon - rc_lon) * sin(rc_lat)  #E_lon_lat
    lhs[2,3] = EARTH_RAD_MM * (rlat - rc_lat)                #E_lat_lat
    
    lhs, rhs
end


"""
    get_block_strain_rates
returns results in nstrain/yr
"""    
function get_block_strain_rates(block_id, vel_df, block_df; min_stations=3,
        ve=:ve, vn=:vn, lon=:lon, lat=:lat)
    vels = vel_df[(vel_df.mov .== block_id), :]
    
    block_lat = block_df[(block_df.fid .== block_id), :lat][1]
    block_lon = block_df[(block_df.fid .== block_id), :lon][1]
    
    if size(vels,1) < min_stations
        strains = (ee=0., en=0., nn=0.)
    else
        lhss = []
        rhss = []
        
        for row in eachrow(vels)
            lhs, rhs = make_strain_rate_equations(
                row[ve], row[vn], row[lon], row[lat], block_lon, block_lat)
            push!(lhss, lhs)
            push!(rhss, rhs)
        end

        LHS = vcat(lhss...)
        RHS = vcat(rhss...)
        
        soln = LHS \ RHS
        
        soln *= 1e9
        
        strains = (ee=soln[1], en=soln[2], nn=soln[3])
    end
    strains
end


function make_proj_matrix(a)

    [cos(a) sin(a); -sin(a) cos(a)]

end


function get_principal_orientation(E)

    if (E[1,1] == 0.) & (E[2,2] == 0.) & (E[1,2] == 0.)
        return 0.
    end

    tan_2theta = (2 * E[2,1]) / (E[1,1] - E[2,2])

    theta = atan(tan_2theta) / 2.
end


function get_principal_stresses(E)
    princ_ang = get_principal_orientation(E)

    P = make_proj_matrix(princ_ang)

    E_prime = P * E * P'

    eig1 = E_prime[1,1]
    eig2 = E_prime[2,2]

    if eig1 > eig2
        eigs = [eig1, eig2]
    else
        eigs = [eig2, eig1]
    end
    eigs
end

end
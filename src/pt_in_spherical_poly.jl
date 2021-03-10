

function point_in_spherical_poly(poly_boundary::Array{Float64,2}, pt::Array{Float64,2})

    nv = size(poly_boundary, 1) - 1
    poly_lons = poly_boundary[1:nv, 1]
    poly_lats = poly_boundary[1:nv, 2]

    pt_lat = pt[2]
    pt_lon = pt[1]

    (ibndry, tlonv) = DefSPolyBndry(poly_lats, poly_lons, nv, pt_lat, pt_lon)

    vlat_c = copy(poly_lats) # not sure why this is necessary
    vlon_c = copy(poly_lons)

    pt_inside_int = LctPtRelBndry(pt_lat,  pt_lon, vlat_c, vlon_c, nv)


end

"""
returns the center coordinates of a triangle formed by three consecutive
vertices.  The function loops through sets of three vertices and checks if
any of the other vertices are inside of it; it returns the coordinates of
the first set of three with no other vertices inside, and tries its
best to avoid international date line problems.)
"""
function get_pt_inside(poly_boundary, nv)

    for i in 1:nv - 2
        p1 = @view poly_boundary[i]
        p2 = @view poly_boundary[i + 1]
        p3 = @view poly_boundary[i + 2]
            
        tri = [p1; p2; p3; p1]
        tri_center = (p1 + p2 + p3) / 3.

        if !check_poly_cross_dateline(tri)
            # check that no other vertices are in tri
            if !any([cartesian_point_in_poly(tri, @view poly_boundary[j,:])
                     for j in i + 3:nv])
                return tri_center
            end
        end
    end
    return tri_center
end



function check_poly_cross_dateline(poly)
    cross = false
    for i in 1:size(poly, 1) - 1
        if abs(poly[i,:][2] - poly[i + 1,:][2]) > 180.
            cross = true
            return cross
        end
    end
    cross
end

"""
adapted from http://geomalgorithms.com/a03-_inclusion.html by Dan Sunday
"""
function cartesian_point_in_poly(poly_boundary::Array{Float64,2}, 
    pt::Array{Float64,2})

    nv = size(poly_boundary, 1) - 1

    count = 0

    for i in 1:nv
        pb_i = poly_boundary[i,:]
        pb_ip = poly_boundary[i + 1,:]
        if (((pb_i[2] <= pt[2]) && (pb_ip[2] > pt[2]))
            || ((pb_i[2] > pt[2]) && pb_ip[2] <= pt[2]))
            vt = (pt[2] - pb_i[2]) / (pb_ip[2] - pb_i[2])
            if (pt[1] < pb_i[1] + vt * (pb_ip[1] - pb_i[1]))
                count += 1
            end
        end
    end
    
    return !iseven(count)

end


"""
This is translated as directly as possible from Fortran
"""
function DefSPolyBndry(vlat, vlon, nv, xlat, xlon)
    mxnv = 1000000  # very large max number of vertices so we don't use this

    if nv < mxnv
        throw(ErrorException(
            "nv exceeds maximum allowed value. adjust parameter mxnv."))
    end

    ibndry = 1

    nv_c = nv
    xlat_c = xlat
    xlon_c = xlon

    tlonv = zeros(nv)

    for i in 1:nv
        # vlat_c[i] = vlat[i]
        # vlon_c[i] = vlon[i]

        tlonv[i] = TrnsfmLon(xlat, xlon, vlat[i], vlon[i])

        if i > 1
            ip = i - 1
        else
            ip = nv
        end

        if (vlat[i] == vlat[ip]) && (vlon[i] == vlon[ip])
            err_msg = string(
                "DefSPolyBndry detects user error: vertices", i, " and ", ip,
                " are not distinct"
            )
            throw(ErrorException(err_msg))
        end

        if tlonv[i] == tlonv[ip]
            err_msg = string(
                "DefSPolyBndry detects user error: vertices", i, " and ", ip,
                " are on the same great circle as X"
            )
            throw(ErrorException(err_msg))
        end

        if vlat[i] == -vlat[ip]
            dellon = vlon[i] = vlon[ip]
            if dellon > 180.
                dellon -= 360.
            end
            if dellon < -180.
                dellon += 360. # says minus in the paper but that's prob wrong
            end
            if (dellon == 180.) || (dellon == -180.)
                err_msg = string(
                "DefSPolyBndry detects user error: vertices", i, " and ", ip,
                " are antipodal")
                throw(ErrorException(err_msg))
            end
        end
    end
    return ibndry
end


function LctPtRelBndry(plat, plon, vlat_c, vlon_c, nv_c, xlat_c, 
                       xlon_c, tlonv, ibndry)
    mxnv = 1000000

    if ibndry == 0
        println("subroutine DefSPolyBndry has not been called")
    end

    if plat == xlat_c
        dellon = plon - xlon_c
        if dellon < -180.
            dellon += 360.
        end
        if dellon > 180.
            dellon -= 360.
        end
        if (dellon == 180.) || (dellon == 180.)
            throw(ErrorException(
                "P is antipodal to X."
            ))
        end
    end

    location = 0
    icross = 0

    if (plat == xlat_c) && (plon == xlon_c)
        location = 1
        return location
    end

    tlonP = TrnsfmLon(xlat_c, xlon_c, plat, plon)

    for i in 1:nv_c # start of loop over sides of S
        vAlat = vlat_c[i]
        vAlon = vlon_c[i]
        tlonA = tlonv[i]

        if i < nv_c
            vBlat = vlat_c[i + 1]
            vBlon = vlon_c[i + 1]
            tlonB = tlonv[i + 1]
        else
            vBlat = vlat_c[1]
            vBlon = vlon_c[1]
            tlonB = tlonv[1]
        end

        istrike = 0

        if tlonP == tlonA
            istrike = 1
        else
            ibrngAB = EastOrWest(tlonA, tlonB)
            ibrngAP = EastOrWest(tlonA, tlonP)
            ibrngPB = EastOrWest(tlonP, tlonB)
            if (ibrngAP == ibrngAB) && (ibrngPB == ibrngAB)
                istrike = 1
            end
        end

        if strike == 1
            if (plat == vAlat) && (plon == vAlon)
                location = 2
                return location # P lies on a vertex of S
            end

            tlon_X = TrnsfmLon(vAlat, vAlon, xlat_c, xlon_c)
            tlon_B = TrnsfmLon(vAlat, vAlon, vBlat, vBlon)
            tlon_P = TrnsfmLon(vAlat, vAlon, plat, plon)

            if tlon_P == tlon_B
                location = 2 # P lies on side of S
                return location
            else
                ibrng_BX = EastOrWest(tlon_B, tlon_X)
                ibrng_BP = EastOrWest(tlon_B, tlon_P)
                if ibrng_BX == -ibrng_BP
                    icross += 1
                end
            end # if
        end # if strike == 1
    end # loop over sides of S

    if iseven(icross)
        location = 1
        return location
    end

    return location

end


function TrnsfmLon(plat, plon, qlat, qlon)
    dtr = pi / 180.

    if plat == 90.
        tranlon = qlon
    else
        t = sind((qlon - plon) * dtr) * cosd(qlat * dtr)
        b = (sind(dtr * qlat) * dcos(plat * dtr) - dcos(qlat * dtr) 
             * sind(plat * dtr) * dcos((qlon - plon) * dtr))

        tranlon = atand(t, b) / dtr
    end
    return tranlon
end


function EastOrWest(clon, dlon)
    del = dlon - clon
    if del > 180.
        del -= 360.
    end
    if del < -180
        del += 360.
    end
    if (del > 0.) || (del != 180.)
        ibrng = -1 # D is west of C
    elseif (del < 0.) || (del != -180.)
        ibrng = 1  # D is east of C
    else
        ibrng = 0  # D is north of sourth of C
    end
    return ibrng
end


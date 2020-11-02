using Test
    
using PyPlot


using Oiler

function test_tri_merc_1()

    tri = Oiler.Tri([81.36, 28.57, -0.],
                    [83.18, 29.05, -25.0],
                    [83.68, 27.56, -0.]
                    )
    
    lon_1 = (81.36 +83.18 + 83.68) / 3.
    lat_1 = (28.57 + 29.05 + 27.56) / 3.

    xp, yp = Oiler.Tris.tri_merc(tri, [lon_1], [lat_1])

    @test isapprox(xp, [135018.41711149618, 0.0, 178067.76749487148, 
                        226987.48383961705])
    @test isapprox(yp, [2.8813557431270294e6, 2.9009189229907184e6,
                        2.9542403022861904e6, 2.7895155119443485e6])
end


function test_tri_strike_slip()

    tri = Oiler.Tri([0., 0., -10.],
                    [1.5, 1.5, -1.],
                    [0., 3., -10.])

    site_lon = [0.35]
    site_lat = [1.5]

    xp, yp = Oiler.Tris.tri_merc(tri, site_lon, site_lat)
    xg, yg = xp[1:1], yp[1:1] # gnss

    zg = zeros(length(xg))

    tx1, ty1 = xp[end - 2], yp[end - 2]
    tx2, ty2 = xp[end - 1], yp[end - 1]
    tx3, ty3 = xp[end], yp[end]

    tz1, tz2, tz3 = tri.p1[3] * 1000., tri.p2[3] * 1000., tri.p3[3] * 1000.

    tp1 = [tx1 ty1 tz1]
    tp2 = [tx2 ty2 tz2]
    tp3 = [tx3 ty3 tz3]


    ss_disp = 10. #1.0e-3 # mm / yr in m
    ds_disp = 0.0 # mm / yr in m
    ts_disp = 0. # no opening or closing of fault

    ue, un, uv = Oiler.TD.TDdispHS(xg, yg, zg, tp1, tp2, tp3, ss_disp, ds_disp)
    uv = zeros(size(uv)) # only considering horizontal velocities

    println(ue, " ", un)

    figure()

    plot([tx1 tx2 tx3 tx1], [ty1 ty2 ty3 ty1], "ro-")
    quiver(xg, yg, ue, un)

    figure()
    plot([tri.p1[1] tri.p2[1] tri.p3[1] tri.p1[1]],
         [tri.p1[2] tri.p2[2] tri.p3[2] tri.p1[2]], "bo-")

    quiver(site_lon, site_lat, ue, un)

    show()

end

#test_tri_dip_slip()



function test_calc_tri_effects_single_tri_1()
    
    tri = Oiler.Tri([0., 0., -10.],
                    [1.5, 1.5, -1.],
                    [0., 3., -10.])

    site_lons = [0.35, 0.25, -0.5]
    site_lats = [1.5, 1.1, 2.3]

    Oiler.Elastic.arrange_tri_partials(
        Oiler.Elastic.calc_tri_effects_single_tri(tri, site_lons, site_lats)...)

end


function test_calc_tri_effects_1()
    tri1 = Oiler.Tri([0., 0., -10.],
                     [1.5, 1.5, -1.],
                     [0., 3., -10.])

    tri2 = Oiler.Tri([0., 3., -10.],
                     [1.5, 1.5, -1.],
                     [1.5, 4.5, -1.])

    site_lons = [0.35, 0.25, -0.5]
    site_lats = [1.5, 1.1, 2.3]

    tris = [tri1, tri2]

    hcat(
        collect(
    
    [Oiler.Elastic.arrange_tri_partials(
        Oiler.Elastic.calc_tri_effects_single_tri(
            tri, site_lons, site_lats)...
        )
    for tri in tris])...
    )
end

#@testset "test tris.jl" begin
#    test_tri_merc_1()
#end

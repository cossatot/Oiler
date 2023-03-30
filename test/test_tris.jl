using Test
    
using PyPlot


using Oiler


function test_get_tri_center()
    tri = Oiler.Tris.Tri(;
                    p1=[0., 0., -10.],
                    p2=[1.5, 1.5, -1.],
                    p3=[0., 3., -10.],
    )
    center = Oiler.Tris.get_tri_center(tri)
    tc_answer = [0.5000510776240583, 1.5001218561577172, -7.0]
    @test isapprox(center, tc_answer)
end


function test_tri_merc_1()

    tri = Oiler.Tris.Tri(;
                    p1=[81.36, 28.57, -0.],
                    p2=[83.18, 29.05, -25.0],
                    p3=[83.68, 27.56, -0.],
                    )
    
    lon_1 = (81.36 + 83.18 + 83.68) / 3.
    lat_1 = (28.57 + 29.05 + 27.56) / 3.

    xp, yp = Oiler.Tris.tri_merc(tri, [lon_1], [lat_1])

    # @test isapprox(xp, [135018.41711149618, 0.0, 178067.76749487148, 
    #                    226987.48383961705])
    # @test isapprox(yp, [2.8813557431270294e6, 2.9009189229907184e6,
    #                    2.9542403022861904e6, 2.7895155119443485e6])
    xp_ans = [-4.399787063193057e6, -4.528173558463124e6, -4.343462673448272e6, -4.326825954910386e6]
    yp_ans = [3.436665555981521e6, 3.483122001542546e6, 3.499694142114799e6, 3.3282431615781095e6]

    @test isapprox(xp, xp_ans)
    @test isapprox(yp, yp_ans)
end


function test_get_tri_strike_dip()
    pt1 = [0., 0., 0.]
    pt2 = [0., 1., 0.]
    pt3 = [1., 0.5, -110.]
    tri_1 = Oiler.Tris.Tri(p1=pt1, p2=pt2, p3=pt3)
    strike_1, dip_1 = Oiler.Tris.get_tri_strike_dip(tri_1)
    @test strike_1 == 0.
    @test isapprox(dip_1, 44.65803022775554)

    pt4 = [0.5, 0.5, 0.]
    pt5 = [1., 0., -110.]
    
    tri_2 = Oiler.Tris.Tri(p1=pt1, p2=pt4, p3=pt5)
    strike_2, dip_2 = Oiler.Tris.get_tri_strike_dip(tri_2)
    @test isapprox(strike_2, 45.19132459503624)
    @test isapprox(dip_2, 54.50372439999114)
end


function test_tri_strike_slip()

    # Use PyPlot for this

    tri = Oiler.Tris.Tri(;
                    p1=[0., 0., -10.],
                    p2=[1.5, 1.5, -1.],
                    p3=[0., 3., -10.],
    )

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


    ss_disp = 10. # 1.0e-3 # mm / yr in m
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

# test_tri_dip_slip()



function test_calc_tri_effects_single_tri_1()
    
    tri = Oiler.Tris.Tri(p1=[0., 0., -10.],
                    p2=[1.5, 1.5, -1.],
                    p3=[0., 3., -10.])

    site_lons = [0.35, 0.25, -0.5]
    site_lats = [1.5, 1.1, 2.3]

    Oiler.Elastic.arrange_tri_partials(
        Oiler.Elastic.calc_tri_effects_single_tri(tri, site_lons, site_lats)...)

end


function test_point_ordering()

    p1 = [0.1, 0.2, -10.]
    p2 = [1.5, 1.5, -1.]
    p3 = [0.1, 3., -10.]

    tri1 = Oiler.Tris.Tri(p1=p1, p2=p2, p3=p3)
    tri2 = Oiler.Tris.Tri(p1=p2, p2=p3, p3=p1)
    tri3 = Oiler.Tris.Tri(p1=p3, p2=p1, p3=p2)
    tri4 = Oiler.Tris.Tri(p1=p3, p2=p2, p3=p1)
    tri5 = Oiler.Tris.Tri(p1=p1, p2=p3, p3=p2)
    tri6 = Oiler.Tris.Tri(p1=p2, p2=p1, p3=p3)
    
    site_lons = [0.35, 0.25, -0.5]
    site_lats = [1.5, 1.1, 2.3]

    outs = [Oiler.Elastic.arrange_tri_partials(
        Oiler.Elastic.calc_tri_effects_single_tri(tri, site_lons, site_lats)...)
        for tri in (tri1, tri2, tri3, tri4, tri5, tri6)]
            
end


function test_calc_tri_effects_1()
    tri1 = Oiler.Tri(p1=[0., 0., -10.],
                     p2=[1.5, 1.5, -1.],
                     p3=[0., 3., -10.])

    tri2 = Oiler.Tri(p1=[0., 3., -10.],
                     p2=[1.5, 1.5, -1.],
                     p3=[1.5, 4.5, -1.])

    site_lons = [0.35, 0.25, -0.5]
    site_lats = [1.5, 1.1, 2.3]

    tris = [tri1, tri2]

    hcat(
        collect(
            [Oiler.Elastic.arrange_tri_partials(
                Oiler.Elastic.calc_tri_effects_single_tri(
                    tri, site_lons, site_lats)...
                )
            for tri in tris]
        )...
    )
end


function test_get_tri_strike_line_1()
    p1 = [0., 1., -1.]
    p2 = [1., 0.5, -0.5]
    p3 = [0., 0., 0.]

    strike_line = Oiler.Tris.get_tri_strike_line(p1, p2, p3)
    @test isapprox(strike_line[1], [1.0, 0.5, -0.5])
    @test isapprox(strike_line[2], [0.0, 0.5000000000000001, -0.5])
end



function test_centroid_distance()
    tri1 = Oiler.Tri(p1=[0., 0., -10.],
                     p2=[1.5, 1.5, -1.],
                     p3=[0., 3., -10.],
                     name="t1")

    tri2 = Oiler.Tri(p1=[0., 3., -10.],
                     p2=[1.5, 1.5, -1.],
                     p3=[1.5, 4.5, -1.], name="t2")

    cd = Oiler.Tris.tri_centroid_distance(tri1, tri2)
    @test isapprox(cd, 175.67806244060736)
end


function test_check_tri_adjacence_1()
    tri1 = Oiler.Tri(p1=[0., 0., -10.],
                     p2=[1.5, 1.5, -1.],
                     p3=[0., 3., -10.],
                     name="t1")

    tri2 = Oiler.Tri(p1=[0., 3., -10.],
                     p2=[1.5, 1.5, -1.],
                     p3=[1.5, 4.5, -1.], name="t2")

    @test Oiler.Tris.check_tri_adjacence(tri1, tri2) == true
end

function test_get_tri_adjacence_dict()
    p1 = [0., 0., 0.]
    p2 = [1.5, 1.5, -10.]
    p3 = [0., 3., 0.]
    p4 = [1.5, 5.5, -10.]
    p5 = [0., 7., 0.]
    p6 = [1.5, 7., -10.]
    p7 = [3., 7., -20.]
    p8 = [3., 3., -20.]
    p9 = [3., 0., -20.]
    p10 = [1.5, 0., -10.]

    tri1 =  Oiler.Tri(p1=p1, p2=p2, p3=p3, name="t1")
    tri2 =  Oiler.Tri(p1=p2, p2=p3, p3=p4, name="t2")
    tri3 =  Oiler.Tri(p1=p3, p2=p4, p3=p5, name="t3")
    tri4 =  Oiler.Tri(p1=p4, p2=p5, p3=p6, name="t4")
    tri5 =  Oiler.Tri(p1=p4, p2=p6, p3=p7, name="t5")
    tri6 =  Oiler.Tri(p1=p4, p2=p7, p3=p8, name="t6")
    tri7 =  Oiler.Tri(p1=p4, p2=p2, p3=p8, name="t7")
    tri8 =  Oiler.Tri(p1=p2, p2=p8, p3=p9, name="t8")
    tri9 =  Oiler.Tri(p1=p10, p2=p2, p3=p9, name="t9")
    tri10 = Oiler.Tri(p1=p1, p2=p2, p3=p10, name="t10")

    tris = [tri1, tri2, tri3, tri4, tri5, tri6, tri7, tri8, tri9, tri10]

    tri_adj_dict_2_pts_no_self_adj = Dict(
        "t3"  => ["t2", "t4"],
        "t7"  => ["t2", "t6", "t8"],
        "t6"  => ["t5", "t7"],
        "t5"  => ["t4", "t6"],
        "t9"  => ["t8", "t10"],
        "t10" => ["t1", "t9"],
        "t2"  => ["t1", "t3", "t7"],
        "t1"  => ["t2", "t10"],
        "t4"  => ["t3", "t5"],
        "t8"  => ["t7", "t9"]
    )

    tri_adj_dict_2_pts_self_adj = Dict(
        "t3"  => ["t2","t3" ,  "t4"],
        "t7"  => [ "t2", "t6", "t7" ,"t8"],
        "t6"  => [ "t5", "t6" ,"t7"],
        "t5"  => [ "t4", "t5" ,"t6"],
        "t9"  => [ "t8","t9" , "t10"],
        "t10" => ["t1", "t9","t10"],
        "t2"  => [ "t1","t2" , "t3", "t7"],
        "t1"  => ["t1" , "t2", "t10"],
        "t4"  => ["t3", "t4" , "t5"],
        "t8"  => [ "t7","t8" , "t9"]
    )

    tri_adj_dict_1_pt_no_self_adj = Dict(
        "t3"  => ["t1", "t2", "t4", "t5", "t6", "t7"],
        "t7"  => ["t1", "t2", "t3", "t4", "t5", "t6", "t8", "t9", "t10"],
        "t6"  => ["t2", "t3", "t4", "t5", "t7", "t8"],
        "t5"  => ["t2", "t3", "t4", "t6", "t7"],
        "t9"  => ["t1", "t2", "t7", "t8", "t10"],
        "t10" => ["t1", "t2", "t7", "t8", "t9"],
        "t2"  => ["t1", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10"],
        "t1"  => ["t2", "t3", "t7", "t8", "t9", "t10"],
        "t4"  => ["t2", "t3", "t5", "t6", "t7"],
        "t8"  => ["t1", "t2", "t6", "t7", "t9", "t10"]
    )

    tri_adj_dict_1_pt_self_adj = Dict(
        "t3"  => ["t1", "t2", "t3" ,"t4", "t5", "t6", "t7"],
        "t7"  => ["t1", "t2", "t3", "t4", "t5", "t6", "t7" ,"t8", "t9", "t10"],
        "t6"  => ["t2", "t3", "t4", "t5", "t6" ,"t7", "t8"],
        "t5"  => ["t2", "t3", "t4","t5" , "t6", "t7"],
        "t9"  => ["t1", "t2", "t7", "t8", "t9" ,"t10"],
        "t10" => ["t1", "t2", "t7", "t8", "t9", "t10"],
        "t2"  => ["t1", "t2" ,"t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10"],
        "t1"  => ["t1" ,"t2", "t3", "t7", "t8", "t9", "t10"],
        "t4"  => ["t2", "t3","t4" , "t5", "t6", "t7"],
        "t8"  => ["t1", "t2", "t6", "t7", "t8" ,"t9", "t10"]
    )

    tri_adj_dict_2_pts_no_self_adj_ = Oiler.Tris.get_tri_adjacence_dict(tris;
        n_common_pts=2, self_adjacence=false)
    
    tri_adj_dict_2_pts_self_adj_ = Oiler.Tris.get_tri_adjacence_dict(tris;
        n_common_pts=2, self_adjacence=true)
    
    tri_adj_dict_1_pt_no_self_adj_ = Oiler.Tris.get_tri_adjacence_dict(tris;
        n_common_pts=1, self_adjacence=false)
    
    tri_adj_dict_1_pt_self_adj_ = Oiler.Tris.get_tri_adjacence_dict(tris;
        n_common_pts=1, self_adjacence=true)
    
    @test tri_adj_dict_1_pt_no_self_adj == tri_adj_dict_1_pt_no_self_adj_
    @test tri_adj_dict_1_pt_self_adj == tri_adj_dict_1_pt_self_adj_
    @test tri_adj_dict_2_pts_no_self_adj == tri_adj_dict_2_pts_no_self_adj_
    @test tri_adj_dict_2_pts_self_adj == tri_adj_dict_2_pts_self_adj_
    
end

@testset "test tris.jl" begin
    test_tri_merc_1()
    test_get_tri_strike_dip()
    test_get_tri_strike_line_1()
    test_centroid_distance()
    test_get_tri_adjacence_dict()
    test_check_tri_adjacence_1()
end

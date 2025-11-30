using Test

import Oiler


function test_gc_dist()
    lax_lon = -(118.0 + 24.0 / 60.0)
    lax_lat = 33.0 + 57.0 / 60.0

    jfk_lon = -(73.0 + 47.0 / 60.0)
    jfk_lat = 40.0 + 38.0 / 60.0

    @test Oiler.gc_distance(lax_lon, lax_lat, jfk_lon, jfk_lat) == 3972.857776250374
end


function test_azimuth_1()
    lax_lon = -(118.0 + 24.0 / 60.0)
    lax_lat = 33.0 + 57.0 / 60.0

    jfk_lon = -(73.0 + 47.0 / 60.0)
    jfk_lat = 40.0 + 38.0 / 60.0

    @test Oiler.azimuth(lax_lon, lax_lat, jfk_lon, jfk_lat) == 65.89216655274531
end


function test_average_azimuth()
    quien_sabe_coords = [-121.37575 36.93768
        -121.33534 36.89357
        -121.29735 36.85596
        -121.25795 36.81047
        -121.2073 36.75499]

    @test Oiler.average_azimuth(quien_sabe_coords[:, 1], quien_sabe_coords[:, 2]) ==
          143.56497333710072
end


function test_beta_from_dip_rake()
    rakes = collect(-180.:45.:180.)
    #rakes[1] = -179.999

    flat_betas = [180.0, 135.0, 90.0, 45.0, 0.0, -45.0, -90.0, -135.0, 180.0] 

    @test flat_betas == Oiler.Geom.beta_from_dip_rake.(0., rakes)

    betas_45 = [180.0 
        144.73561031724535
        90.0
        35.26438968275466
        0.0
        -35.26438968275466
        -90.0
        -144.73561031724535
        180.0]

    @test isapprox(betas_45, Oiler.Geom.beta_from_dip_rake.(45, rakes))

end


function test_az_from_strike_dip_rake()
    rakes = collect(-180.:45.:180.)
    
    azs_0_0 = [180., 135., 90., 45., 0., 315., 270., 225., 180.]

   @test isapprox(azs_0_0, Oiler.Geom.az_from_strike_dip_rake.(0., 0.0001, rakes))

    azs_0_45 = [180.0, 144.73561031724535, 90.0, 35.264389682754654, 0.0, 
                324.73561031724535, 270.0, 215.26438968275465, 180.0]

    @test isapprox(azs_0_45, Oiler.Geom.az_from_strike_dip_rake.(0., 45., rakes))

    azs_45_45 = Oiler.Geom.wrap360.(azs_0_45 .+ 45.0)

    @test isapprox(azs_45_45, Oiler.Geom.az_from_strike_dip_rake.(45., 45., rakes))

end


function test_rake_from_az_strike_dip()
    rakes = collect((-179:180))
    #azs = collect(0.:359.)
    strikes = collect(0.:45.:270.)
    dips = [1, 20, 45, 60, 89]

    for rake in rakes
        for strike in strikes
            for dip in dips
                az = Oiler.Geom.az_from_strike_dip_rake(strike, dip, rake)
                test_rake = Oiler.Geom.rake_from_az_strike_dip(az, strike, dip)
                test_pass = isapprox(rake, test_rake)
                if !test_pass
                    beta = az-strike
                    println(
                    "rake: $rake, test_rake: $test_rake, az: $az, strike: $strike, beta: $beta")
                end
                @test isapprox(rake, test_rake)
            end
        end
    end
end
                    

function test_oblique_merc_2()
    lat1 = 0.00001
    lat2 = -0.00001
    lon1 = -1.0
    lon2 = 1.0

    lons = [0.25; 0.35]
    lats = [0.25; 0.25]

    x, y = Oiler.oblique_merc(lons, lats, lon1, lat1, lon2, lat2)

end

function test_terminal_coords_from_bearing_dist_1()
    @test Oiler.Geom.terminal_coords_from_bearing_dist(0.0, 0.0, 0.0, 10.0) ==
          (0.0, 0.08993216059187306)
end

function test_terminal_coords_from_bearing_dist_2()
    @test Oiler.Geom.terminal_coords_from_bearing_dist(0.0, 89.5, 0.0, 100.0) ==
          (180.0, 89.60067839408168)
end

function test_terminal_coords_from_bearing_dist_3()
    @test Oiler.Geom.terminal_coords_from_bearing_dist(0.0, 0.0, 90.0, 10.0) ==
          (0.08993216059187306, 0.0)
end


function test_polyline_length()
    trace = [0.0 0.0; 0.0 0.5; 0.0 1.0]
    @test Oiler.Geom.polyline_length(trace) == 111.19492664455873
end


function test_point_line_distance_1()
    point = [1.0, 1.0]
    line_start = [0.0, 0.0]
    line_end = [2.0, 0.0]
    dist = Oiler.Geom.point_line_distance(point, line_start, line_end)

    @test dist == 1.0
end


function test_point_line_distance_2()
    point = [0.0, 0.0]
    line_start = [0.0, 0.0]
    line_end = [2.0, 0.0]
    dist = Oiler.Geom.point_line_distance(point, line_start, line_end)

    @test dist == 0.0
end


function test_simplify_polyline_vert()
    polyline = [0.0 0.0; 0.0 1.0; 0.0 2.0; 0.0 3.0; 0.0 4.0]
    simp_line = [0.0 0.0; 0.0 4.0]
    test_simp_line = Oiler.Geom.simplify_polyline(polyline, 0.0)

    @test test_simp_line == simp_line
end


function test_simplify_polyline_rosettacode()
    polyline = [0.0 0.0
        1.0 0.1
        2.0 -0.1
        3.0 5.0
        4.0 6.0
        5.0 7.0
        6.0 8.1
        7.0 9.0
        8.0 9.0
        9.0 9.0]
    simp_line = [0.0 0.0; 2.0 -0.1; 3.0 5.0; 7.0 9.0; 9.0 9.0]
    test_simp_line = Oiler.Geom.simplify_polyline(polyline, 1.0)
    @test test_simp_line == simp_line
end


function test_break_polyline_equal_1()
    tr = [0.0 0.0; 3.0 3.0; 9.0 9.0; 10.0 10.0]
    @test Oiler.Geom.break_polyline_equal(tr, 1) == [tr]
end

function test_break_polyline_equal_2()
    tr = [0.0 0.0; 3.0 3.0; 9.0 9.0; 10.0 10.0]
    @test Oiler.Geom.break_polyline_equal(tr, 2) ==
          [[0.0 0.0; 3.0 3.0; 4.980085962376654 5.000771545105226],
        [4.980085962376654 5.000771545105226; 9.0 9.0; 10.0 10.0]]

end

function test_break_polyline_equal_3()
    tr = [0.0 0.0; 3.0 3.0; 9.0 9.0; 10.0 10.0]
    @test Oiler.Geom.break_polyline_equal(tr, 4) ==
          [[0.0 0.0; 2.493526173498305 2.4945825672649673],
        [2.493526173498305 2.4945825672649673; 3.0 3.0; 4.980085962376654 5.000771545105226],
        [4.980085962376654 5.000771545105226; 7.481540598820497 7.50169718205662],
        [7.481540598820497 7.50169718205662; 9.0 9.0; 10.0 10.0]]
end


function test_check_winding_order_array_2_1()
    aa = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.0 0.0]
    @test Oiler.Geom.check_winding_order(aa) == 1
end

function test_check_winding_order_array_arrays_1()
    bb = [[0.0 0.0], [1.0 0.0], [1.0 1.0], [0.0 1.0], [0.0 0.0]]
    @test Oiler.Geom.check_winding_order(bb) == 1
end


function test_strike_dip_from_3_pts()
    pt1 = [0.0, 0.0, 0.0]
    pt2 = [0.0, 1.0, 0.0]
    pt3 = [1.0, 0.5, -1.0]

    strike_1, dip_1 = Oiler.Geom.strike_dip_from_3_pts(pt1, pt2, pt3)
    @test strike_1 == 0.0
    @test dip_1 == 45.0

    pt4 = [0.5, 0.5, 0.0]
    pt5 = [1.0, 0.0, -1.0]

    strike_2, dip_2 = Oiler.Geom.strike_dip_from_3_pts(pt1, pt4, pt5)
    @test strike_2 == 45.0
    @test isapprox(dip_2, 54.735610317245346)

    pt6 = [-1.0, -1.0, 0.0]
    pt7 = [0.0, -1.0, -1.0]

    strike_3, dip_3 = Oiler.Geom.strike_dip_from_3_pts(pt1, pt6, pt7)
    @test strike_3 == 45.0
    @test isapprox(dip_3, 54.735610317245346)

    pt8 = [0.0, -1.0, 1.0]
    strike_4, dip_4 = Oiler.Geom.strike_dip_from_3_pts(pt1, pt6, pt8)
    @test strike_4 == 225.0
    @test isapprox(dip_4, 54.735610317245346)

    pt9 = [0.0, 0.5, -1.0]
    # fails due to SingularException.  Could add a fix?
    # strike_5, dip_5 = Oiler.Geom.strike_dip_from_3_pts(pt1, pt2, pt9)
    # @test strike_5 == 0.
    # @test dip_5 == 90.

end



@testset "test geom.jl" begin
    test_gc_dist()
    test_azimuth_1()
    test_average_azimuth()
    test_beta_from_dip_rake()
    test_az_from_strike_dip_rake()
    test_rake_from_az_strike_dip()
    test_oblique_merc_2()
    test_terminal_coords_from_bearing_dist_1()
    test_terminal_coords_from_bearing_dist_2()
    test_terminal_coords_from_bearing_dist_3()
    test_polyline_length()
    test_point_line_distance_1()
    test_point_line_distance_2()
    test_simplify_polyline_vert()
    test_simplify_polyline_rosettacode()
    test_break_polyline_equal_1()
    test_break_polyline_equal_2()
    test_break_polyline_equal_3()
    test_strike_dip_from_3_pts()

end


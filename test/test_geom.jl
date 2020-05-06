using Test

import Oiler

lax_lon = 118. + 24. / 60.;
lax_lat = 33. + 57. / 60.;

jfk_lat = 40. + 38. / 60.;
jfk_lon = 73. + 47. / 60.;


quien_sabe_coords = [-121.37575 36.93768;
                     -121.33534 36.89357;
                     -121.29735 36.85596;
                     -121.25795 36.81047;
                     -121.2073 36.75499];

@test Oiler.gc_distance(lax_lon, lax_lat, jfk_lon, jfk_lat) == 3972.857776250374

#@test Oiler.azimuth(lax_lon, lax_lat, jfk_lon, jfk_lat) == 77.11787999624288

#@test Oiler.average_azimuth(quien_sabe_coords[:,1], quien_sabe_coords[:,2]) ==
#   323.8414207049782

function test_oblique_merc_2()
    lat1 = 0.00001
    lat2 = -0.00001
    lon1 = -1.
    lon2 = 1.

    lons = [0.25; 0.35]
    lats = [0.25; 0.25]

    x, y = Oiler.oblique_merc(lons, lats, lon1, lat1, lon2, lat2)

end

test_oblique_merc_2()


function test_point_line_distance_1()
    point = [1., 1.]
    line_start = [0., 0.]
    line_end = [2., 0.]
    dist = Oiler.Geom.point_line_distance(point, line_start, line_end)

    @test dist == 1.0
end

test_point_line_distance_1()

function test_point_line_distance_2()
    point = [0., 0.]
    line_start = [0., 0.]
    line_end = [2., 0.]
    dist = Oiler.Geom.point_line_distance(point, line_start, line_end)

    @test dist == 0.0
end

test_point_line_distance_2()

function test_simplify_polyline_vert()
    polyline = [0. 0.; 0. 1.; 0. 2.; 0. 3.; 0. 4.]
    simp_line = [0. 0.; 0. 4.]
    test_simp_line = Oiler.Geom.simplify_polyline(polyline, 0.)
    @test test_simp_line == simp_line
end

test_simplify_polyline_vert()

function test_simplify_polyline_rosettacode()
    polyline = [0. 0.; 1. 0.1; 2. -0.1; 3. 5.; 4. 6.; 
                5. 7.; 6. 8.1; 7. 9.; 8. 9.; 9. 9.]
    simp_line = [0. 0.; 2. -0.1; 3. 5.; 7. 9.; 9. 9.]
    test_simp_line = Oiler.Geom.simplify_polyline(polyline, 1.)
    @test test_simp_line == simp_line
end

test_simplify_polyline_rosettacode()
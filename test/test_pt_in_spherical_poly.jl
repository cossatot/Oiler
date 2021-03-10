using Test

include("../src/pt_in_spherical_poly.jl")

function test_cartesian_point_in_poly_1()
    sq = [0. 0.; 0. 1.; 1. 1.; 1. 0.; 0. 0.]
    pt = [0.5 0.5]

    @test cartesian_point_in_poly(sq, pt) == true
end

function test_cartesian_point_in_poly_2()
    sq = [0. 0.; 0. 1.; 1. 1.; 1. 0.; 0. 0.]
    pt = [1.5 1.5]

    @test cartesian_point_in_poly(sq, pt) == false
end

function test_poly_cross_dateline_cross_1()
    poly = [1. -179.; 1. 179.; 2. 179.; 1. -179]
    @test check_poly_cross_dateline(poly) == true
end


function test_poly_cross_dateline_no_cross_1()
    poly = [1. 179.; 1. 179.1; 2. 179.; 1. 179]
    @test check_poly_cross_dateline(poly) == false
end

function test_get_pt_inside_poly_in_1()

    poly = [[0. 0.], [1. 0.], [1. 1.], [0. 1.], [0. 0.]]

    println(get_pt_inside_poly(poly))
end



@testset "test pt_in_spherical_poly" begin
   test_cartesian_point_in_poly_1()
   test_cartesian_point_in_poly_2()
   test_poly_cross_dateline_cross_1()
   test_poly_cross_dateline_no_cross_1()
end
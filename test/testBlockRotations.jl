include("../src/BlockRotations.jl")

using Test

@test build_PvGb_Deg(0., 0.) == [ 0.0  -6.371e6  0.0;
                                  0.0   0.0      6.371e6;
                                  0.0   0.0      0.0]

@test build_Pv_deg(0., 0.) == [ -0.0  -0.0   1.0;
                                -0.0   1.0   0.0;
                                -1.0  -0.0  -0.0]

@test build_Gb_deg(0., 0.) == [  0.0   0.0      -0.0;
                                -0.0   0.0       6.371e6;
                                 0.0  -6.371e6   0.0]

function test_euler_pole_conversion_1()
    lond = 36.
    latd = 50.
    rotr = 0.2

    pole_s = EulerPoleSphere(lond, latd, rotr)
    pole_c = euler_pole_sphere_to_cart(pole_s)

    @test euler_pole_cart_to_sphere(pole_c) == pole_s
end

@test test_euler_pole_conversion_1()

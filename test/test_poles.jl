import Oiler

using Test

function test_euler_pole_conversion_1()
    lon = 36.
    lat = 50.
    rotr = 0.2

    pole_s = Oiler.RotationPoles.PoleSphere(lon = lon, lat = lat, rotrate = rotr)
    pole_c = Oiler.RotationPoles.pole_sphere_to_cart(pole_s)

    pole_sc = Oiler.RotationPoles.pole_cart_to_sphere(pole_c) 
    @test isapprox(pole_sc.lon, pole_s.lon)
    @test isapprox(pole_sc.lat, pole_s.lat)
    @test isapprox(pole_sc.rotrate, pole_s.rotrate)
end

@testset "test_euler_pole_conversion_1" begin
    test_euler_pole_conversion_1()
end
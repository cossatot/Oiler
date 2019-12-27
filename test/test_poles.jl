import Oiler

using Test

function test_euler_pole_conversion_1()
    lond = 36.
    latd = 50.
    rotr = 0.2

    pole_s = Oiler.RotationPoles.PoleSphere(lond=lond, latd=latd, rotrate=rotr)
    pole_c = Oiler.RotationPoles.pole_sphere_to_cart(pole_s)

    pole_sc = Oiler.RotationPoles.pole_cart_to_sphere(pole_c) 
    @test isapprox(pole_sc.lond, pole_s.lond)
    @test isapprox(pole_sc.latd, pole_s.latd)
    @test isapprox(pole_sc.rotrate, pole_s.rotrate)
end

test_euler_pole_conversion_1()
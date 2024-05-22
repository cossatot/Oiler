import Oiler

using Test

function test_euler_pole_conversion_1()
    lon = 36.
    lat = 50.
    rotr = 0.2

    pole_s = Oiler.RotationPoles.PoleSphere(lon=lon, lat=lat, rotrate=rotr)
    pole_c = Oiler.RotationPoles.pole_sphere_to_cart(pole_s)

    pole_sc = Oiler.RotationPoles.pole_cart_to_sphere(pole_c) 
    @test isapprox(pole_sc.lon, pole_s.lon)
    @test isapprox(pole_sc.lat, pole_s.lat)
    @test isapprox(pole_sc.rotrate, pole_s.rotrate)
end


function test_minus()
    p1 = Oiler.PoleCart(x=100., y=200., z=-300., fix="me", mov="it")
    minus_p1 = -p1

    @test minus_p1.x == -100.0
    @test minus_p1.y == -200.0
    @test minus_p1.z == 300.0
    @test minus_p1.ex == 0.0
    @test minus_p1.ey == 0.0
    @test minus_p1.ez == 0.0
    @test minus_p1.fix == "it"
    @test minus_p1.mov == "me"
end


function test_err_cart_to_sphere()

    pole = Oiler.PoleCart(
        x=4.942785033008153e-10,
        y=3.0887336517757523e-9,
        z=-2.864680725578047e-9,
        ex=1.4821768706642342e-10,
        ey=1.4821747112934177e-10,
        ez=1.482185081352683e-10,
        fix="ca",
        mov="na"
    )

    elon, elat, erotrate = Oiler.RotationPoles.err_cart_to_sphere(
        pole.x, pole.y, pole.z, pole.ex, pole.ey, pole.ez, pole.cxy,
        pole.cxz, pole.cyz)

    @test isapprox(elon, 0.15338052234690125)
    @test isapprox(elat, 0.06996314752420013)
    @test isapprox(erotrate, 1.2587072976554877e-12)
end


function test_cart_to_sphere_w_errors()
    pole = Oiler.PoleCart(
        x=4.942785033008153e-10,
        y=3.0887336517757523e-9,
        z=-2.864680725578047e-9,
        ex=1.4821768706642342e-10,
        ey=1.4821747112934177e-10,
        ez=1.482185081352683e-10,
        fix="ca",
        mov="na"
    )

    pole_sphere = Oiler.RotationPoles.pole_cart_to_sphere(pole)

    @test isapprox(pole_sphere.lon, 80.90825589356842)
    @test isapprox(pole_sphere.lat, -42.483737633837634)
    @test isapprox(pole_sphere.rotrate, 0.24302450801360484)
    @test isapprox(pole_sphere.elon, 0.15338052234690125)
    @test isapprox(pole_sphere.elat, 0.06996314752420013)
    @test isapprox(pole_sphere.erotrate, 1.2587072976554877e-12)
    @test pole_sphere.fix == "ca"
    @test pole_sphere.mov == "na"
end


function test_err_cart_to_sphere_and_back()
    pole_cart = Oiler.PoleCart(
        x=1.7595239620687005e-9,
        y=-7.509137368732926e-9,
        z=-4.11009580217663e-9,
        ex=1.4887106312185033e-9,
        ey=9.481715572889885e-9,
        ez=6.380368177374067e-9,
        cxy=1.405005292532337e-17,
        cxz=9.457727139873025e-18,
        cyz=6.047061530260166e-17,
        fix="1111",
        mov="6"
    )

 
    pole_sphere = Oiler.pole_cart_to_sphere(pole_cart)
    pole_cart_back = Oiler.pole_sphere_to_cart(pole_sphere)

    @test isapprox(pole_cart.x,   pole_cart_back.x)
    @test isapprox(pole_cart.y,   pole_cart_back.y)
    @test isapprox(pole_cart.z,   pole_cart_back.z)
    @test isapprox(pole_cart.ex,  pole_cart_back.ex)
    @test isapprox(pole_cart.ey,  pole_cart_back.ey)
    @test isapprox(pole_cart.ez,  pole_cart_back.ez)
    @test pole_cart.fix == pole_cart_back.fix
    @test pole_cart.mov == pole_cart_back.mov
end


@testset "test_euler_pole_conversion_1" begin
    test_euler_pole_conversion_1()
    test_err_cart_to_sphere()
    test_cart_to_sphere_w_errors()
    test_minus()
    # test_err_cart_to_sphere_and_back()
end
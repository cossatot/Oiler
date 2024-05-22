import Oiler

using Test

function test_predict_block_vel_1()
    PvGb = Oiler.BlockRotations.build_PvGb_deg(-13.98, -52.17)
    
    pole = Oiler.PoleSphere(lat=9.3, lon=-41.7, rotrate=0.15, 
        mov="an", fix="af");
    pole_cart = Oiler.pole_sphere_to_cart(pole)

    V = PvGb * [pole_cart.x; pole_cart.y; pole_cart.z]

    @test isapprox(V, 
        [13.16176174592891; 7.656387837884508; -2.299548549083007e-16])
end


function test_direct_solve_1_vec()
    PvGb = Oiler.BlockRotations.build_PvGb_deg(-13.98, -52.17)
    V = [13.16176174592891; 7.656387837884508; -2.299548549083007e-16]
    p = PvGb \ V
    pole = Oiler.pole_cart_to_sphere(Oiler.PoleCart(x=p[1], y=p[2], z=p[3]))

    @test isapprox(pole.lon, -41.7)
    @test isapprox(pole.lat, 9.3)
    @test isapprox(pole.rotrate, 0.15)
end


function test_predict_block_vel_2()
    
    ca_na_pole = Oiler.PoleSphere(lat=64.9, lon=-110.5, rotrate=0.214, 
               elon=10.0, elat=4.0, erotrate=0.03, mov="ca", fix="na")

    site_lon = -81.3
    site_lat = 12.53

    Oiler.BlockRotations.predict_block_vel(site_lon, site_lat, pole)
   
end



function test_make_block_PvGb_from_vels_1_vel()
end

@testset "test block_rotations.jl" begin
    test_predict_block_vel_1()
    test_direct_solve_1_vec()
end

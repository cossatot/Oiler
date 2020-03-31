import Oiler

using Test

function test_predict_block_vel_1()
    PvGb = Oiler.BlockRotations.build_PvGb_deg(-13.98, -52.17)
    
    pole = Oiler.PoleSphere(latd=9.3, lond=-41.7, rotrate=0.15, 
        mov="an", fix="af");
    pole_cart = Oiler.pole_sphere_to_cart(pole)

    V = PvGb * [pole_cart.x; pole_cart.y; pole_cart.z]

    @test isapprox(V, 
        [13.16176174592891; 7.656387837884508; -2.299548549083007e-16])
end

test_predict_block_vel_1()

function test_direct_solve_1_vec()
    PvGb = Oiler.BlockRotations.build_PvGb_deg(-13.98, -52.17)
    V = [13.16176174592891; 7.656387837884508; -2.299548549083007e-16]
    p = PvGb \ V
    pole = Oiler.pole_cart_to_sphere(Oiler.PoleCart(x=p[1], y=p[2], z=p[3]))

    @test isapprox(pole.lond, -41.7)
    @test isapprox(pole.latd, 9.3)
    @test isapprox(pole.rotrate, 0.15)
end

test_direct_solve_1_vec()


function test_make_block_PvGb_from_vels_1_vel()
end


function test_solve_block_invs_from_vel_groups_1_vel()
    V = Oiler.VelocityVectorSphere(lond=-13.98, latd=-52.17,
        ve = 13.16176174592891, vn = 7.656387837884508, fix="af", mov="an")

    vg = Oiler.group_vels_by_fix_mov([V])

    poles = Oiler.solve_block_invs_from_vel_groups(vg)
end

#test_solve_block_invs_from_vel_groups_1_vel()
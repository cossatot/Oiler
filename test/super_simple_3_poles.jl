using Revise

using Test

using Oiler

na_pa = Oiler.PoleSphere(lat = 49.3002664, lon = -76.0123145, 
    rotrate = 0.790821495, mov = "na", fix = "pa");
jf_pa = Oiler.PoleSphere(lat = -0.6, lon = 37.8, 
    rotrate = 0.625, mov = "jf", fix = "pa");
na_jf = na_pa - jf_pa

p_na_1 = (-121.8, 46.8)
p_na_2 = (-122.6, 45.5)

p_jf_1 = (-127.4, 45.3)
p_jf_2 = (-126.0, 41.3)

p_ca_1 = (-117.1, 44.0)
p_ca_2 = (-126.0, 47.2)


na_1, na_2 = Oiler.predict_block_vels([p_na_1[1], p_na_2[1]], 
                                      [p_na_1[2], p_na_2[2]],
                                      na_pa)

jf_1, jf_2 = Oiler.predict_block_vels([p_jf_1[1], p_jf_2[1]],
                                      [p_jf_1[2], p_jf_2[2]],
                                      jf_pa)

ca_1, ca_2 = Oiler.predict_block_vels([p_ca_1[1], p_ca_2[1]], 
                                      [p_ca_1[2], p_ca_2[2]],
                                      na_jf)


# modify these to add weights                                      
na_1 = Oiler.VelocityVectorSphere( lon = -121.8, lat = 46.8, ve = 16.48835031187386,
    vn = -41.100626589504365, vu = 0.0, ee = 1.0, en = 1.0, eu = 1.0,
    fix = "pa", mov = "na", name = "na_1")

na_2 = Oiler.VelocityVectorSphere( lon = -122.6, lat = 45.5, ve = 18.61975328133655,
    vn = -41.65490927154829, vu = 0.0, ee = 1.0, en = 1.0, eu = 1.0,
    fix = "pa", mov = "na", name = "na_2")

jf_1 = Oiler.VelocityVectorSphere( lon = -127.4, lat = 45.3, ve = 47.24492282732255,
    vn = -17.751696806912953, vu = 0.0, ee = 1.0, en = 1.0, eu = 1.0,
    fix = "pa", mov = "jf", name = "jf_1")

jf_2 = Oiler.VelocityVectorSphere( lon = -126.0, lat = 41.3, ve = 43.497621077223826,
    vn = -19.387934120552465, vu = 0.0, ee = 1.0, en = 1.0, eu = 1.0,
    fix = "pa", mov = "jf", name = "jf_2")

ca_1 = Oiler.VelocityVectorSphere( lon = -117.1, lat = 44.0, ve = -25.258056124365822,
    vn = -8.207145561346188, vu = 0.0, ee = 1.0, en = 1.0, eu = 1.0,
    fix = "jf", mov = "na", name = "ca_1")

ca_2 = Oiler.VelocityVectorSphere( lon = -126.0, lat = 47.2, ve = -30.2250503313903,
    vn = -24.53080399521363, vu = 0.0, ee = 1.0, en = 1.0, eu = 1.0,
    fix = "jf", mov = "na", name = "ca_2")


all_vels = [na_1, na_2, jf_1, jf_2, ca_1, ca_2]

vel_groups = Oiler.group_vels_by_fix_mov(all_vels)
no_weight_poles = Oiler.solve_block_invs_from_vel_groups(vel_groups, 
                                                         weighted = false)
nwp = Dict(k => Oiler.pole_cart_to_sphere(v) for (k, v) in no_weight_poles)

weight_poles = Oiler.solve_block_invs_from_vel_groups(vel_groups, 
                                                         weighted = true)
wp = Dict(k => Oiler.pole_cart_to_sphere(v) for (k, v) in weight_poles)

function test_poles_equal()
    for (k, p) in wp
        @test isapprox(p.lon, nwp[k].lon; rtol = 1e-5)
        @test isapprox(p.lat, nwp[k].lat; rtol = 1e-5)
        @test isapprox(p.rotrate, nwp[k].rotrate; rtol = 1e-5)
    end
end

@testset "super_simple_3_poles.jl pole equality" begin
    test_poles_equal()
end
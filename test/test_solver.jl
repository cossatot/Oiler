using Test

include("../src/solver.jl")

using Oiler


function test_build_constraint_matrix()
    keys = [("a", "b"), ("b", "c"), ("c", "a"), ("r", "a"), ("r", "c")]
end


function test_build_constraint_matrices()
end


function test_weight_from_error_nonzero()
end


function test_weight_from_error_zero()
end

function test_build_weight_vector_from_vel()
end


function test_build_weight_vector_from_vels_default_zero_weight()
    vv = Oiler.VelocityVectorSphere(lon = 0., lat = 0., ve = 1., vn = 1.,  en = 2.,
    ee = 3.)
    
    @test Oiler.Solver.build_weight_vector_from_vels([vv]) == [2.; 3.; 0.]
end

test_build_weight_vector_from_vels_default_zero_weight()


function test_build_weight_vector_from_vels_equal_zero_weights()
    vv = Oiler.VelocityVectorSphere(lon = 0., lat = 0., ve = 1., vn = 1., ee = 0.,
    en = 0.)
    
    ww = Oiler.Solver.build_weight_vector_from_vels([vv])

    @test ww == [1; 1; 1]
end

test_build_weight_vector_from_vels_equal_zero_weights()


function test_make_block_PvGb_from_vels()

    vels = [
        Oiler.VelocityVectorSphere(lon = -13.98, lat = -52.17, ve = 1., vn = 1.,
        ee = 0., en = 0.)]

end


function test_solve_block_invs_from_vel_groups_1_vel()
    V = Oiler.VelocityVectorSphere(lon = -13.98, lat = -52.17,
        ve = 13.16176174592891, vn = 7.656387837884508, fix = "af", mov = "an")

    vg = Oiler.group_vels_by_fix_mov([V])

    poles = Oiler.solve_block_invs_from_vel_groups(vg)

    
end

# test_solve_block_invs_from_vel_groups_1_vel()



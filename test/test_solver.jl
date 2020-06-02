using Test

using SparseArrays

using Oiler





function test_build_constraint_matrix()
    keys = [("a", "b"), ("b", "c"), ("c", "a"), ("r", "a"), ("r", "c")]
    cycles = Oiler.Utils.find_vel_cycles(keys)
    cycle = cycles[1]
    cm = Oiler.Solver.build_constraint_matrix(cycle, keys)
    cm_answer = Oiler.Utils.dict_to_sparse(
        Dict(
                    (1, 7) => 1.0,
                    (2, 8) => 1.0,
                    (3, 9) => 1.0,
                    (1, 10) => -1.0,
                    (2, 11) => -1.0,
                    (3, 12) => -1.0,
                    (1, 13) => 1.0,
                    (2, 14) => 1.0,
                    (3, 15) => 1.0
                    ))
    @test cm == cm_answer
end


function test_build_constraint_matrices()
    keys = [("a", "b"), ("b", "c"), ("c", "a"), ("r", "a"), ("r", "c")]
    cycles = Oiler.Utils.find_vel_cycles(keys)
    cm = Oiler.Solver.build_constraint_matrices(cycles, keys)
    cm_answer = Oiler.Utils.dict_to_sparse(
        Dict(
    (1,  1)  => 1.0,
    (2,  2)  => 1.0,
    (3,  3)  => 1.0,
    (1,  4)  => 1.0,
    (2,  5)  => 1.0,
    (3,  6)  => 1.0,
    (1,  7)  => 1.0,
    (4,  7)  => 1.0,
    (2,  8)  => 1.0,
    (5,  8)  => 1.0,
    (3,  9)  => 1.0,
    (6,  9)  => 1.0,
    (4, 10)  => -1.0,
    (5, 11)  => -1.0,
    (6, 12)  => -1.0,
    (4, 13)  => 1.0,
    (5, 14)  => 1.0,
    (6, 15)  => 1.0
        )
    )

    @test cm == cm_answer
end


function test_weight_from_error_nonzero()
    @test Oiler.Solver.weight_from_error(3.) == 1 / 9.
end


function test_weight_from_error_zero()
    @test Oiler.Solver.weight_from_error(0.) == 1.0e20
end


function test_weight_from_error_zero_specified()
    @test Oiler.Solver.weight_from_error(0.; zero_err_weight = 2.) == 0.25
end


function test_build_weight_vector_from_vel()
    vv = Oiler.VelocityVectorSphere(lon = 0., lat = 0., ve = 1., vn = 1., 
    en = 2., ee = 3.)
    
    @test Oiler.Solver.build_weight_vector_from_vel(vv) == [1 / 9.; 0.25; 1.0e20]

end


function test_build_weight_vector_from_vels_default_zero_weight()
    vv = Oiler.VelocityVectorSphere(lon = 0., lat = 0., ve = 1., vn = 1., 
        en = 2., ee = 3.)
    ww = Oiler.VelocityVectorSphere(lon = 0., lat = 0., ve = 1., vn = 1., 
        en = 1., ee = 1.)

    weight = Oiler.Solver.build_weight_vector_from_vels([vv, ww])
    weight_answer = [1 / 9.; 0.25; 1.0e20; 1.; 1.; 1.0e20]

    @test weight == weight_answer
end



function test_build_weight_vector_from_vels_equal_zero_weights()
    vv = Oiler.VelocityVectorSphere(lon = 0., lat = 0., ve = 1., vn = 1., ee = 0.,
    en = 0.)
    
    ww = Oiler.Solver.build_weight_vector_from_vels([vv])

    @test ww == [1.; 1.; 1.]
end


function test_build_weight_vectors()
    vv = Oiler.VelocityVectorSphere(lon = 0., lat = 0., ve = 1., vn = 1., 
    en = 2., ee = 3.)
    ww = Oiler.VelocityVectorSphere(lon = 0., lat = 0., ve = 1., vn = 1., 
    en = 1., ee = 1.)

    vel_groups = Dict(("a", "b") => [vv, ww])
    weight = Oiler.Solver.build_weight_vectors(vel_groups)
    weight_answer = Dict(("a", "b") => [1 / 9.; 0.25; 1.0e20; 1.; 1.; 1.0e20])
    @test weight == weight_answer
end



function test_make_block_PvGb_from_vel()

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
@testset "basic Solver tests" begin
    test_build_constraint_matrix()
    test_build_constraint_matrices()
    test_weight_from_error_nonzero()
    test_weight_from_error_zero()
    test_weight_from_error_zero_specified()
    test_build_weight_vector_from_vel()
    test_build_weight_vector_from_vels_default_zero_weight()
    test_build_weight_vector_from_vels_equal_zero_weights()
    test_build_weight_vectors()
end
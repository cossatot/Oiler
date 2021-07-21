using Test

using SparseArrays
using LinearAlgebra

using Oiler

faults_1 = [Oiler.Fault(trace=[0. 0.1; 1. 1.], dip=89., dip_dir="SW", hw="a",
                    fw="c", name="ff1"),
          Oiler.Fault(trace=[5. 5.; 1. 1.], dip=89., dip_dir="SW", hw="a",
                    fw="b", name="ff2"),
          Oiler.Fault(trace=[5. 5.; 1. 1.], dip=89., dip_dir="NE", hw="c",
                    fw="b", name="ff3"),
          Oiler.Fault(trace=[5. 5.; 1. 1.], dip=89., dip_dir="SW", hw="b",
                    fw="c", name="ff4"),
]

vel_groups_1 = Dict(
    # ("a", "b") => [
    # VelocityVectorSphere(lon=0.0,
    #                      lat=0.05,
    #                      ve=0.5,
    #                      vn=0.,
    #                      fix="a",
    #                      mov="b",
    #                      name="f1",
    #                      vel_type="fault"),
    # VelocityVectorSphere(lon=0.05,
    #                      lat=0.,
    #                      ve=0.53,
    #                      vn=0.,
    #                      fix="a",
    #                      mov="b",
    #                      name="f2",
    #                      vel_type="fault")
    # ],

    ("a", "c") => [
    VelocityVectorSphere(lon=0.7,
                         lat=0.05,
                         ve=1.5,
                         vn=0.,
                         fix="c",
                         mov="b",
                         name="f3",
                         vel_type="fault")
],
    ("r", "a") => [
    VelocityVectorSphere(lon=-1.,
                         lat=0.05,
                         ve=0.0,
                         vn=0.,
                         fix="r",
                         mov="a",
                         name="g3",
                         vel_type="GNSS"),

    VelocityVectorSphere(lon=0.05,
                         lat=0.,
                         ve=0.15,
                         vn=0.,
                         fix="r",
                         mov="a",
                         name="g4",
                         vel_type="GNSS")
    ],
    ("r", "c") => [
    VelocityVectorSphere(lon=0.77,
                         lat=0.05,
                         ve=1.9,
                         vn=0.,
                         fix="r",
                         mov="c",
                         name="g5",
                         vel_type="GNSS")
])

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
    @test Oiler.Solver.weight_from_error(0.) == 1.0e-4
end


function test_weight_from_error_zero_specified()
    @test Oiler.Solver.weight_from_error(0.; zero_err_weight=2.) == 0.25
end


function test_build_weight_vector_from_vel()
    vv = Oiler.VelocityVectorSphere(lon=0., lat=0., ve=1., vn=1., 
    en=2., ee=3.)
    
    @test Oiler.Solver.build_weight_vector_from_vel(vv) == [1 / 9.; 0.25; 1.0e-4]

end


function test_build_weight_vector_from_vels_default_zero_weight()
    vv = Oiler.VelocityVectorSphere(lon=0., lat=0., ve=1., vn=1., 
        en=2., ee=3.)
    ww = Oiler.VelocityVectorSphere(lon=0., lat=0., ve=1., vn=1., 
        en=1., ee=1.)

    weight = Oiler.Solver.build_weight_vector_from_vels([vv, ww])
    weight_answer = [1 / 9.; 0.25; 1.0e-4; 1.; 1.; 1.0e-4]

    @test weight == weight_answer
end



function test_build_weight_vector_from_vels_equal_zero_weights()
    vv = Oiler.VelocityVectorSphere(lon=0., lat=0., ve=1., vn=1., ee=0.,
    en=0.)
    ww = Oiler.VelocityVectorSphere(lon=1., lat=0., ve=1., vn=1., ee=0.,
    en=0.)
    
    vel_groups = Dict(("a", "b") => [vv, ww])
    block_matrices = Oiler.Solver.make_block_PvGb_from_vels(vel_groups)

    @test block_matrices["weights"] == [1.; 1.; 1.; 1.; 1.; 1]
end


function test_build_weight_vectors()
    vv = Oiler.VelocityVectorSphere(lon=0., lat=0., ve=1., vn=1., 
    en=2., ee=3.)
    ww = Oiler.VelocityVectorSphere(lon=0., lat=0., ve=1., vn=1., 
    en=1., ee=1.)

    vel_groups = Dict(("a", "b") => [vv, ww])
    weight = Oiler.Solver.build_weight_vectors(vel_groups)
    weight_answer = Dict(("a", "b") => [1 / 9.; 0.25; 1.0e-4; 1.; 1.; 1.0e-4])
    @test weight == weight_answer
end


function test_make_block_PvGb_from_vel()

    vels = [
        Oiler.VelocityVectorSphere(lon=-13.98, lat=-52.17, ve=1., vn=1.,
        ee=0., en=0., fix="ca", mov="na")]

    vel_groups = Oiler.group_vels_by_fix_mov(vels)

    Oiler.Solver.make_block_PvGb_from_vels(vel_groups)

    pvgb = Oiler.Solver.make_block_PvGb_from_vels(vel_groups)

    pvgb_ans = Dict("keys" => [("ca", "na")],
                    "weights" => [1.0, 1.0, 1.0],
                    "PvGb" => sparse([4.88298e9   -1.21565e9  3.90747e9;
                    -1.53913e9   -6.18229e9  0.0;
                    -1.19209e-7   0.0        0.0]))

    @test sort(collect(keys(pvgb))) == sort(collect(keys(pvgb_ans)))
    @test pvgb["weights"] == pvgb_ans["weights"]
    @test isapprox(pvgb["PvGb"], pvgb_ans["PvGb"]; rtol=0.01)
end

function test_make_block_PvGb_from_vels()

    pvgb = Oiler.Solver.make_block_PvGb_from_vels(vel_groups_1)

    pvgb_ans = Dict("keys" => [("a", "c"), ("r", "a"), ("r", "c")],
                    "weights" => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                    "PvGb" => sparse(
    [-5.55933e6 -67923.4 6.371e9 0.0 0.0 0.0 0.0 0.0 0.0;
    7.78345e7 -6.37052e9 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    0.0 0.0 0.0 -5.5589e6 97030.9 6.371e9 0.0 0.0 0.0;
    0.0 0.0 0.0 -1.11189e8 -6.37003e9 1.45519e-11 0.0 0.0 0.0;
    0.0 0.0 0.0 1.45519e-11 0.0 -1.49012e-8 0.0 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 6.371e9 0.0 0.0 0.0;
    0.0 0.0 0.0 5.55975e6 -6.371e9 0.0 0.0 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 -9.31323e-10 0.0 0.0 0.0;
    0.0 0.0 0.0 0.0 0.0 0.0 -5.55924e6 -74715.4 6.371e9;
    0.0 0.0 0.0 0.0 0.0 0.0 8.56175e7 -6.37042e9 0.0;
    0.0 0.0 0.0 0.0 0.0 0.0 -1.45519e-11 0.0 0.0])
    )

    @test sort(collect(keys(pvgb))) == sort(collect(keys(pvgb_ans)))
    @test pvgb["weights"] == pvgb_ans["weights"]
    @test isapprox(pvgb["PvGb"], pvgb_ans["PvGb"]; rtol=0.01)

end


function test_make_block_inversion_matrices_from_vels()
end


function test_make_block_inv_rhs()
end

function test_make_block_inv_lhs_constraints()
end

function test_add_fault_locking_to_PvGb()
end


function test_make_tri_regularization_matrix()
end

function test_make_tri_prior_matrices()
end


function test_add_tris_to_PvGb()

end

function test_weight_inv_matrices()

end


function test_add_equality_constraints_kkt()
end


function test_solve_block_invs_from_vel_groups_1_vel()
    V = Oiler.VelocityVectorSphere(lon=-13.98, lat=-52.17,
        ve=13.16176174592891, vn=7.656387837884508, fix="af", mov="an")

    vg = Oiler.group_vels_by_fix_mov([V])

    poles = Oiler.solve_block_invs_from_vel_groups(vg)

    
end

function test_set_up_block_inv_w_constraints_1()
    block_inv_setup = Oiler.Solver.set_up_block_inv_w_constraints(vel_groups_1)

end


function test_set_up_block_inv_w_constraints_1_tri()
    tri = Oiler.Tri(p1=[-1.01, -1.0, 0.],
                    p2=[-1.51, 2.0, -25.0],
                    p3=[-1.01, 2., 0.])

    block_inv_setup = Oiler.Solver.set_up_block_inv_w_constraints(vel_groups_1,
        tris=[tri])

end

"""
This function contains a set of tests that compare different linear
least-squares solution strategies: ordinary LLS (OLS), weighted LLS (WLS),
equality-constrained LLS (CLS), equality-constrained weighted LLS (CWLS).

The tests ensure that the weights are being handled correctly, that the
constraints are met, and that non-weighted and weighted results converge
when the weights all equal 1.
""" 
function test_solver_strategies()

    # set up
    x_obs = [0.,  1.,    2.,   3.,    4.]
    y_obs = [0.3, 2.005, 2.98, 5.2,   5.]
    y_e =   [1.4, 0.1,   5.,   0.1,   0.1]
    y_w =   map(Oiler.Solver.weight_from_error, y_e)
    y_w1 = ones(size(y_w))

    cm = [1. -1.]
    cm0 = [0. 0.]
    c = 0.

    PvGb = [0. 1.; 1. 1.; 2. 1.; 3. 1.; 4. 1.]

    # no weights
    m, b = PvGb \ y_obs

    # weighted (real weights)
    w_lhs, w_rhs = Oiler.Solver.weight_inv_matrices(PvGb, y_obs, y_w)
    m_w, b_w = w_lhs \ w_rhs
    @test isapprox(m_w, 1.086484313142308)
    @test isapprox(b_w, 1.1695147339514742)
    
    # weighted (equal weights, should be identical to OLS)
    w1_lhs, w1_rhs = Oiler.Solver.weight_inv_matrices(PvGb, y_obs, y_w1)
    m_w1, b_w1 = w1_lhs \ w1_rhs
    @test isapprox(m_w1, m)
    @test isapprox(b_w1, b)

    # constrained least squares
    c_lhs, c_rhs = Oiler.Solver.add_equality_constraints_kkt(PvGb, y_obs, cm)
    m_c, b_c = c_lhs \ c_rhs
    @test isapprox(m_c, b_c)

    # constrained, weighted least squares (real weights)
    cw_lhs, cw_rhs = Oiler.Solver.make_weighted_constrained_lls_matrices(PvGb,
        y_obs, cm, y_w)
    # m_cw, b_cw, _, _, _, _, _, _ = cw_lhs \ cw_rhs
    _, _, _, _, _, _, m_cw, b_cw = cw_lhs \ cw_rhs
    @test isapprox(m_cw, b_cw)
    @test isapprox(m_cw, m_w, atol=0.1)
    @test isapprox(b_cw, b_w, atol=0.1)
    
    ## constrained, weighted least squares (zero constraints, should match WLS)
    ## fails w/ singular value exception
    # cw0_lhs, cw0_rhs = Oiler.Solver.make_weighted_constrained_lls_matrices(PvGb,
    #    y_obs, cm0, y_w)
    # _, _, _, _, _, _, m_cw0, b_cw0 = cw0_lhs \ cw0_rhs
    # @test isapprox(m_w, m_cw0)
    # @test isapprox(b_w, b_cw0)

    # KKT constrained, weighted least squares (equal weights, should match CLS)
    cw1_lhs, cw1_rhs = Oiler.Solver.make_weighted_constrained_kkt_lls_matrices(PvGb,
        y_obs, cm, y_w1)
    m_cw1, b_cw1, _, _, _, _, _, _ = cw1_lhs \ cw1_rhs
    # KKT constrained, weighted least squares (equal weights, should match CLS)
    @test isapprox(m_cw1, b_cw1)
    @test isapprox(m_c, m_cw1)
    @test isapprox(b_c, b_cw1)
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
    test_make_block_PvGb_from_vel()
    test_make_block_PvGb_from_vels()
    # test_add_tris_to_PvGb()
    # test_solve_block_invs_from_vel_groups_1_vel()
    # test_set_up_block_inv_w_constraints_1()
    # test_set_up_block_inv_w_constraints_1_tri()

    test_solver_strategies() 
end
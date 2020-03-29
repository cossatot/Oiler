using Test

include("../src/solver.jl")

using Solver
using Solver: VelocityVectorSphere

function test_build_weight_vector_from_vels()
    vv = Solver.VelocityVectorSphere(lond = 0., latd = 0., ve = 1., vn = 1.,  en = 2.,
    ee = 3.)
    
    @test Solver.build_weight_vector_from_vels([vv]) == [0.5; 0.3333333333333333;
    100000.0]
end


function test_make_block_PvGb_from_vels()

    vels = [VelocityVectorSphere()]



    test_build_weight_vector_from_vels()
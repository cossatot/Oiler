using Test

include("../src/solver.jl")

using Oiler

function test_build_weight_vector_from_vels_default_zero_weight()
    vv = Oiler.VelocityVectorSphere(lond = 0., latd = 0., ve = 1., vn = 1.,  en = 2.,
    ee = 3.)
    
    @test Oiler.Solver.build_weight_vector_from_vels([vv]) == [2.; 3.; 0.]
end

test_build_weight_vector_from_vels_default_zero_weight()


function test_build_weight_vector_from_vels_equal_zero_weights()
    vv = Oiler.VelocityVectorSphere(lond=0., latd=0., ve=1., vn=1., ee=0.,
    en=0.)
    
    ww = Oiler.Solver.build_weight_vector_from_vels([vv])

    @test ww == [1; 1; 1]
end

test_build_weight_vector_from_vels_equal_zero_weights()

function test_make_block_PvGb_from_vels()

    vels = [
        Oiler.VelocityVectorSphere(lond=-13.98, latd=-52.17, ve=1., vn=1.,
        ee=0., en=0.)]

    
end



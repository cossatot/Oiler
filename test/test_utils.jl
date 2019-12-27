using Test

include("../src/utils.jl")

function test_build_weight_vector_from_vels()
    vv = Utils.VelocityVectorSphere(lond=0., latd=0., ve=1., vn=1.,  en=2.,
    ee=3.)
    
    @test Utils.build_weight_vector_from_vels([vv]) == [0.5; 0.3333333333333333;
    100000.0]
end

test_build_weight_vector_from_vels()
using Test

using Oiler

@testset "Oiler.jl unit tests" begin
    include("test_poles.jl")
    include("test_solver.jl")
end

@testset "Oiler.jl integration tests" begin
end
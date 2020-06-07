using Test

using Oiler

@testset "Oiler.jl unit tests" begin
    # work up the heirarchy dependencies
    include("test_poles.jl")
    # include("test_velocities.jl")
    include("test_geom.jl")
    include("test_faults.jl")
    include("test_block_rotations.jl")
    include("test_io.jl")
    include("test_utils.jl")
    include("test_okada.jl")
    # include("test_elastic.jl")  # need to fix many tests
    include("test_solver.jl")
end

@testset "Oiler.jl integration tests" begin
    include("super_simple_3_poles.jl")
end
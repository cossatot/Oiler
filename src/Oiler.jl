module Oiler

include("constants.jl")
using .Constants: EARTH_RAD_KM, EARTH_RAD_MM
export EARTH_RAD_KM, EARTH_RAD_MM

include("poles.jl")
using .RotationPoles:  PoleCart, PoleSphere, add_poles, add_poles, subtract_poles,
    pole_cart_to_sphere, pole_sphere_to_cart
export PoleCart, PoleSphere, add_poles, add_poles, subtract_poles,
    pole_cart_to_sphere, pole_sphere_to_cart

include("velocities.jl")
using .Velocities: VelocityVectorSphere
export VelocityVectorSphere

include("geom.jl")
using .Geom: azimuth, gc_distance, average_azimuth, az_to_angle, angle_to_az,
    angle_difference
export azimuth, gc_distance, average_azimuth, az_to_angle, angle_to_az,
    angle_difference

include("faults.jl")
using .Faults: Fault, fault_to_vel, fault_slip_rate_to_ve_vn
export Fault, fault_to_vel, fault_slip_rate_to_ve_vn

include("block_rotations.jl")
using .BlockRotations: build_PvGb_from_vels, build_vel_column_from_vels, predict_block_vels
export build_PvGb_from_vels, build_vel_column_from_vels, predict_block_vels

include("io.jl")
using .IO: vel_from_row, load_vels_from_csv, group_vels_by_fix_mov
export vel_from_row, load_vels_from_csv, group_vels_by_fix_mov

include("utils.jl")
using .Utils: make_block_PvGb_from_vels, solve_block_invs_from_vel_groups,
    predict_vels_from_poles, solve_for_block_poles_iterative
export make_block_PvGb_from_vels, solve_block_invs_from_vel_groups,
    predict_vels_from_poles, solve_for_block_poles_iterative

end  # module
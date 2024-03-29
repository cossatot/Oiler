module Oiler

include("constants.jl")
using .Constants: EARTH_RAD_KM, EARTH_RAD_MM
export EARTH_RAD_KM, EARTH_RAD_MM

include("stats.jl")
using .Stats

include("poles.jl")
using .RotationPoles:  PoleCart, PoleSphere, add_poles, add_poles, subtract_poles,
    pole_cart_to_sphere, pole_sphere_to_cart
export PoleCart, PoleSphere, add_poles, add_poles, subtract_poles,
    pole_cart_to_sphere, pole_sphere_to_cart

include("velocities.jl")
using .Velocities: VelocityVectorSphere, reverse
export VelocityVectorSphere, reverse

include("geom.jl")
using .Geom: azimuth, gc_distance, average_azimuth, az_to_angle, angle_to_az,
    angle_difference, rotate_velocity, rotate_xy_vec, oblique_merc
export azimuth, gc_distance, average_azimuth, az_to_angle, angle_to_az,
    angle_difference, rotate_velocity, rotate_xy_vec, oblique_merc

include("faults.jl")
using .Faults: Fault, fault_to_vel, fault_slip_rate_to_ve_vn, 
    ve_vn_to_fault_slip_rate, fault_oblique_merc, build_strike_rot_matrix,
    build_velocity_projection_matrix, fault_to_vels
export Fault, fault_to_vel, fault_slip_rate_to_ve_vn,
    ve_vn_to_fault_slip_rate, fault_oblique_merc, build_strike_rot_matrix,
    build_velocity_projection_matrix, fault_to_vels

include("block_rotations.jl")
using .BlockRotations: build_PvGb_from_vels, build_vel_column_from_vels, predict_block_vels
export build_PvGb_from_vels, build_vel_column_from_vels, predict_block_vels

include("bounds.jl")
using .Boundaries

include("io.jl")
using .IO: vel_from_row, load_vels_from_csv, group_vels_by_fix_mov
export vel_from_row, load_vels_from_csv, group_vels_by_fix_mov

include("utils.jl")
using .Utils: predict_vels_from_poles, find_vel_cycles, diagonalize_matrices, 
    get_gnss_vels, get_coords_from_vel_array
export predict_vels_from_poles, find_vel_cycles, diagonalize_matrices, 
    get_gnss_vels, get_coords_from_vel_array

include("okada.jl")
using .Okada:okada
export okada

include("tdisphs.jl")
using .TD

include("tris.jl")
using .Tris

include("elastic.jl")
using .Elastic:fault_to_okada
export fault_to_okada

include("strain.jl")
using .Strain

include("results_analysis.jl")
using .ResultsAnalysis

include("solver.jl")
using .Solver: make_block_PvGb_from_vels, solve_block_invs_from_vel_groups,
    solve_for_block_poles_iterative, 
    make_block_inversion_matrices_from_vels
export make_block_PvGb_from_vels, solve_block_invs_from_vel_groups,
    predict_vels_from_poles, solve_for_block_poles_iterative, find_vel_cycles,
    make_block_inversion_matrices_from_vels


include("plots.jl")
using .Plots

include("write_web_viewer.jl")
using .WebViewer

end  # module

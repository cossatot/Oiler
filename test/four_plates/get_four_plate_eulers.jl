using Revise

using DataFrames
using CSV
using PyPlot
using Random

using Oiler

Random.seed!(69)

vels = Oiler.load_vels_from_csv("./data/fault_vels_err.csv");
vel_groups = Oiler.group_vels_by_fix_mov(vels);

poles = Oiler.solve_block_invs_from_vel_groups(vel_groups);

poles_iter = Oiler.solve_for_block_poles_iterative(vel_groups, 5);

println("blah")
using Revise

using DataFrames
using CSV
using PyPlot

using Oiler

vels = Oiler.load_vels_from_csv("./data/fault_vels_err.csv");
vel_groups = Oiler.group_vels_by_fix_mov(vels);

poles = Oiler.solve_block_invs_from_vel_groups(vel_groups);
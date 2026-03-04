using Oiler
using PyPlot

block_df = Oiler.IO.gis_vec_file_to_df("./test_data/anatolia_locking_test_blocks.geojson")

fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
    "./test_data/anatolia_locking_test_faults.geojson",
    block_df=block_df)

vel_df = Oiler.IO.gis_vec_file_to_df("./test_data/anatolia_locking_test_vels.csv")
gnss_vels = Oiler.IO.make_vels_from_gnss_and_blocks(vel_df, block_df;
    ve=:veast_eu, vn=:vnorth_eu, ee=:seast, en=:snorth, name=:names, fix="eu")

# Inversion 1: all faults
non_fault_bounds1 = Oiler.IO.get_non_fault_block_bounds(block_df, faults)
bv1 = vcat(map(b -> Oiler.Boundaries.boundary_to_vels(b), non_fault_bounds1)...)
vg1 = Oiler.group_vels_by_fix_mov(vcat(fault_vels, bv1, gnss_vels))
results1 = Oiler.Solver.solve_block_invs_from_vel_groups(vg1; faults=faults,
    predict_vels=true)

# Inversion 2: without anf008
faults2 = filter(f -> f.fid != "anf008", faults)
fault_vels2 = Oiler.IO.make_vels_from_faults(faults2)
non_fault_bounds2 = Oiler.IO.get_non_fault_block_bounds(block_df, faults2)
bv2 = vcat(map(b -> Oiler.Boundaries.boundary_to_vels(b), non_fault_bounds2)...)
vg2 = Oiler.group_vels_by_fix_mov(vcat(fault_vels2, bv2, gnss_vels))
results2 = Oiler.Solver.solve_block_invs_from_vel_groups(vg2; faults=faults2,
    predict_vels=true)

# Inversion 3: without anf008, with off-fault locking at 50 km
results3 = Oiler.Solver.solve_block_invs_from_vel_groups(vg2; faults=faults2,
    non_fault_bounds=non_fault_bounds2, off_fault_locking_depth=50.0, predict_vels=true)

fig1 = Oiler.Plots.plot_results_map(results1, vg1, faults)
title("all faults")

fig2 = Oiler.Plots.plot_results_map(results2, vg2, faults2)
title("no anf008")

fig3 = Oiler.Plots.plot_results_map(results3, vg2, faults2)
title("no anf008 + locking 50 km")

#show()

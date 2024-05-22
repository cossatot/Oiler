using Revise

using Test

using DataFrames

using Oiler




test_gj_blocks = "./test_data/io_test_blocks.geojson"
test_gj_faults = "./test_data/io_test_faults.geojson"
test_gj_boundary = "./test_data/io_test_boundary.geojson"
test_gj_gnss_vels = "./test_data/io_test_gnss_vels.geojson"
test_gj_geol_slip_rates = "./test_data/io_test_geol_slip_rates.geojson"

pred_vel_pt_file = "./test_data/io_pt_grid_low_res.geojson"

block_df = Oiler.IO.gis_vec_file_to_df(test_gj_blocks)

fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
    test_gj_faults, block_df=block_df,
    subset_in_bounds=true)

gnss_vel_df = Oiler.IO.gis_vec_file_to_df(test_gj_gnss_vels)
gnss_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gnss_vel_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err)


slip_rate_df = Oiler.IO.gis_vec_file_to_df(test_gj_geol_slip_rates)
slip_rate_df, slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(slip_rate_df,
    fault_df)

vels = vcat(fault_vels,
    gnss_vels,
    slip_rate_vels)

vel_groups = Oiler.group_vels_by_fix_mov(vels)

results = Oiler.solve_block_invs_from_vel_groups(vel_groups; faults=faults,
    predict_vels=true, pred_se=true)

#map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults)
#rates_fig = Oiler.Plots.plot_slip_rate_fig(slip_rate_df, slip_rate_vels,
#    fault_df, results)




vel_pt_df = Oiler.IO.load_pred_vels_from_pt_file(pred_vel_pt_file)
pred_vel_pts = Oiler.IO.make_vels_for_pred_from_pt_file(pred_vel_pt_file,
    block_df)


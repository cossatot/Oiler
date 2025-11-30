using Revise

using Oiler

using PyPlot


AVES = Oiler.VelocityVectorSphere(lon = -63.62, lat = 15.67, ve = 13.2, vn =
14.2, ee = 2.9, en = 2.2, fix = "na", mov = "ca");
ROJO = Oiler.VelocityVectorSphere(lon = -71.67, lat = 17.90, ve = 13.5, vn =
6.90, ee = 3.0, en = 2.6, fix = "na", mov = "ca");
SANA = Oiler.VelocityVectorSphere(lon = -81.73, lat = 12.53, ve = 14.8, vn =
7.4, ee = 2.7, en = 2.5, fix = "na", mov = "ca");
CRO1 = Oiler.VelocityVectorSphere(lon = -64.58, lat = 17.76, ve = 9.80, vn =
13.3, ee = 2.0, en = 2.1, fix = "na", mov = "ca");

vels = [AVES; ROJO; SANA; CRO1];

vel_groups = Oiler.group_vels_by_fix_mov(vels);

weighted_stuff = Oiler.Solver.set_up_block_inv_w_constraints(vel_groups, 
    weighted=true) 

nonweighted_stuff = Oiler.Solver.set_up_block_inv_w_constraints(vel_groups, 
    weighted=false) 

wbi = Oiler.Solver.make_block_inversion_matrices_from_vels(vel_groups)
bi = Oiler.Solver.make_block_inversion_matrices_from_vels(vel_groups)

poles = Oiler.solve_block_invs_from_vel_groups(vel_groups; weighted = true);
npoles = Oiler.solve_block_invs_from_vel_groups(vel_groups; weighted = false);

pred_vels = Oiler.predict_block_vels(vel_groups[("na", "ca")], poles[("na", "ca")])
npred_vels = Oiler.predict_block_vels(vel_groups[("na", "ca")], npoles[("na", "ca")])

ve_pred = [pv.ve for pv in pred_vels];
vn_pred = [pv.vn for pv in pred_vels];

nve_pred = [pv.ve for pv in npred_vels];
nvn_pred = [pv.vn for pv in npred_vels];

lons, lats = Oiler.get_coords_from_vel_array(npred_vels)

ve_obs, vn_obs = [v.ve for v in vels], [v.vn for v in vels]

figure()
quiver(lons, lats, nve_pred, nvn_pred, color = "black", scale = 100.)
quiver(lons, lats, ve_pred, vn_pred, color = "blue", scale = 100.)
quiver(lons, lats, ve_obs, vn_obs, color = "red", scale = 100.)
show()
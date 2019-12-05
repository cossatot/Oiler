using Revise

using DataFrames
using CSV
using PyPlot

using Oiler

vels = Oiler.load_vels_from_csv("./data/fault_vels.csv");
vel_groups = Oiler.group_vels_by_fix_mov(vels);

af_eu_vels = vel_groups[("eu","af")]
na_af_vels = vel_groups[("af","na")]
na_eu_vels = vel_groups[("eu","na")]
na_sa_vels = vel_groups[("sa","na")]
af_sa_vels = vel_groups[("sa","af")]

af_eu_PvGb = Oiler.build_PvGb_from_vels(af_eu_vels);
na_af_PvGb = Oiler.build_PvGb_from_vels(na_af_vels);
na_eu_PvGb = Oiler.build_PvGb_from_vels(na_eu_vels);
na_sa_PvGb = Oiler.build_PvGb_from_vels(na_sa_vels);
af_sa_PvGb = Oiler.build_PvGb_from_vels(af_sa_vels);


big_PvGb = Oiler.diagonalize_matrices((af_eu_PvGb, na_af_PvGb, na_eu_PvGb, na_sa_PvGb,
                                 af_sa_PvGb));


big_vels = Oiler.build_vel_column_from_vels([af_eu_vels; na_af_vels; na_eu_vels;
                                             na_sa_vels; af_sa_vels]);

constraint_set_11 = [1 0 0 1 0 0 -1 -0 -0 0 0 0 0 0 0];
constraint_set_12 = [0 1 0 0 1 0 -0 -1 -0 0 0 0 0 0 0];
constraint_set_13 = [0 0 1 0 0 1 -0 -0 -1 0 0 0 0 0 0];

constraint_set_21 = [0 0 0 1 0 0 0 0 0 -1 -0 -0 1 0 0];
constraint_set_22 = [0 0 0 0 1 0 0 0 0 -0 -1 -0 0 1 0];
constraint_set_23 = [0 0 0 0 0 1 0 0 0 -0 -0 -1 0 0 1];

constraint_lhs = vcat(constraint_set_11, constraint_set_12, constraint_set_13,
                      constraint_set_21, constraint_set_22, constraint_set_23);

constraint_lhs = convert(Array{AbstractFloat}, constraint_lhs);

constraint_rhs = zeros(size(constraint_lhs)[1]);

m, n  = size(big_PvGb);
p = size(constraint_rhs)[1]

lhs_term_1 = 2 * big_PvGb' * big_PvGb

lhs = [lhs_term_1 constraint_lhs'; constraint_lhs zeros(p, p)]

rhs = [2 * big_PvGb' * big_vels; constraint_rhs];

kkt_sol = lhs \ rhs
omegas = kkt_sol[1:n]

#omegas = big_PvGb \ big_vels;

af_eu_pole = Oiler.EulerPoleCart(x = omegas[1], y = omegas[2], z = omegas[3],
                                 mov = "af", fix = "eu");
na_af_pole = Oiler.EulerPoleCart(x = omegas[4], y = omegas[5], z = omegas[6],
                                 mov = "na", fix = "af");
na_eu_pole = Oiler.EulerPoleCart(x = omegas[7], y = omegas[8], z = omegas[9],
                                 mov = "na", fix = "eu");
na_sa_pole = Oiler.EulerPoleCart(x = omegas[10], y = omegas[11], z = omegas[12],
                                 mov = "na", fix = "sa");
af_sa_pole = Oiler.EulerPoleCart(x = omegas[13], y = omegas[14], z = omegas[15],
                                 mov = "af", fix = "sa");

sa_eu_pole_1 = Oiler.euler_pole_cart_to_sphere(Oiler.add_poles(na_sa_pole, -na_eu_pole));
sa_eu_pole_2 = Oiler.euler_pole_cart_to_sphere(Oiler.add_poles(na_af_pole, af_eu_pole));

na_eu_lats = [vel.latd for vel in na_eu_vels];
na_eu_lons = [vel.lond for vel in na_eu_vels];

na_eu_pred = Oiler.predict_block_vels(na_eu_lons, na_eu_lats, na_eu_pole);
na_eu_pole_af = Oiler.add_poles(af_eu_pole, na_af_pole)

na_eu_pred_af = Oiler.predict_block_vels(na_eu_lons, na_eu_lats, na_eu_pole_af);


figure()
quiver(na_eu_lons, na_eu_lats, [v.ve for v in na_eu_vels], 
       [v.vn for v in na_eu_vels], color = "blue", scale = 100)
quiver(na_eu_lons, na_eu_lats, [v.ve for v in na_eu_pred], 
       [v.vn for v in na_eu_pred], color = "red", scale = 100)
quiver(na_eu_lons, na_eu_lats, [v.ve for v in na_eu_pred_af], 
       [v.vn for v in na_eu_pred_af], color = "green", scale = 100)
show()
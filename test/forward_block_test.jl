include("../src/BlockRotations.jl")
using .BlockRotations
const br = BlockRotations

using PyPlot


AVES = br.VelocityVectorSph(-63.62, 15.67, 13.2, 14.2, 0., 0., 0., 0.);
ROJO = br.VelocityVectorSph(-71.67, 17.9, 13.5, 6.9, 0., 0., 0., 0.);
SANA = br.VelocityVectorSph(-81.73, 12.53, 14.8, 7.4, 0., 0., 0., 0.);
CRO1 = br.VelocityVectorSph(-64.58, 17.76, 9.8, 13.3, 0., 0., 0., 0.);

vels = [AVES; ROJO; SANA; CRO1];

A = br.build_PvGb_from_vels(vels);

V = br.build_vel_column_from_vels(vels);

omega_hat = A \ V;

euler_cart = br.EulerPoleCart(omega_hat[1], omega_hat[2],
                           omega_hat[3]);

euler_sphere = br.euler_pole_cart_to_sphere(euler_cart)

xs = [AVES.lond ROJO.lond SANA.lond CRO1.lond];
ys = [AVES.latd ROJO.latd SANA.latd CRO1.latd];

pred_vels = br.predict_block_vels(xs, ys, euler_sphere);

ve_obs = [AVES.ve ROJO.ve SANA.ve CRO1.ve];
vn_obs = [AVES.vn ROJO.vn SANA.vn CRO1.vn];

#ve_pred = [(pv_aves * gb_aves * omega_hat)[1]
#           (pv_rojo * gb_rojo * omega_hat)[1]
#           (pv_sana * gb_sana * omega_hat)[1]
#           (pv_cro1 * gb_cro1 * omega_hat)[1]];
#
#vn_pred = [(pv_aves * gb_aves * omega_hat)[2]
#           (pv_rojo * gb_rojo * omega_hat)[2]
#           (pv_sana * gb_sana * omega_hat)[2]
#           (pv_cro1 * gb_cro1 * omega_hat)[2]];
#
#

ve_pred = [pv.ve for pv in pred_vels];
vn_pred = [pv.vn for pv in pred_vels];

figure()
quiver(xs, ys, ve_pred, vn_pred, color = "blue", scale = 100.)
quiver(xs, ys, ve_obs, vn_obs, color = "red", scale = 100.)
show()
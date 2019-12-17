using Revise

using Oiler

using PyPlot


AVES = Oiler.VelocityVectorSphere(lond = -63.62, latd = 15.67, ve = 13.2, vn = 14.2);
ROJO = Oiler.VelocityVectorSphere(lond = -71.67, latd = 17.90, ve = 13.5, vn = 6.90);
SANA = Oiler.VelocityVectorSphere(lond = -81.73, latd = 12.53, ve = 14.8, vn = 7.4);
CRO1 = Oiler.VelocityVectorSphere(lond = -64.58, latd = 17.76, ve = 9.80, vn = 13.3);

vels = [AVES; ROJO; SANA; CRO1];

A = Oiler.build_PvGb_from_vels(vels);

V = Oiler.build_vel_column_from_vels(vels);

omega_hat = A \ V;

euler_cart = Oiler.EulerPoleCart(x = omega_hat[1], y = omega_hat[2], z = omega_hat[3],
                                 fix = "", mov = "");

euler_sphere = Oiler.euler_pole_cart_to_sphere(euler_cart)

xs = [AVES.lond ROJO.lond SANA.lond CRO1.lond];
ys = [AVES.latd ROJO.latd SANA.latd CRO1.latd];

pred_vels = Oiler.predict_block_vels(xs, ys, euler_sphere);

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
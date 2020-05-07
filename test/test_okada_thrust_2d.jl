using Revise
using PyPlot

using Oiler

pole_sph = Oiler.PoleSphere(lon = -10, lat = 10., rotrate = 1. / (19. * sqrt(2)));
pole_cart = Oiler.pole_sphere_to_cart(pole_sph);
pole_vec = [pole_cart.x; pole_cart.y; pole_cart.z];


thrust_trace = [-10. 10.; 10. -10.];
thrust = Oiler.Fault(trace = thrust_trace, dip_dir = "SW", dip = 20., lsd = 20.);

vlons = collect(-1.:0.01:1.);
vlats = collect(-1.:0.01:1.);

f_pvgb = Oiler.BlockRotations.build_PvGb_deg(0., 0.);

partials = Oiler.Elastic.calc_locking_effects_per_fault(thrust, vlons, vlats);

pred_vecs = [part * pole_vec for part in partials];


pe, pn, pv = [v[1] for v in pred_vecs], [v[2] for v in pred_vecs], [v[3] for v in pred_vecs];

p_hor = sqrt.(pe.^2 .+ pn.^2);

dist = sqrt.((-1. .- vlons ).^2 + (-1 .- vlats).^2);

figure()
# plot(dist, p_hor)
plot(vlons, pe, "b")
plot(vlons, pn, "r")
plot(vlons, pv, "g")

figure()
plot([-10.; 10.], [10.; -10], "k")
plot(vlons, vlats, "r.")
# quiver(vlons, vlats, pe, pn)
show()

# f_pvgb * pole_vec


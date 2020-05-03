using Revise

using Oiler

using PyPlot


pole = Oiler.PoleCart(x = -4.000480063406787e-10, y = -3.686194245405466e-9,
                      z = 3.100995264274769e-9, fix = "na", mov = "ca");

# thrust_trace = [-81.905471 12.945880; -81.554797 12.114126];
# thrust = Oiler.Fault(trace = thrust_trace, dip_dir = "SW", dip = 20., hw = "na",
#                      fw = "ca")

strike_slip_trace = [-81.5 13.75; -80.650 14.091; -79.798 14.430];
ss = Oiler.Fault(trace = strike_slip_trace, dip_dir = "S", dip = 89., hw = "ca", fw = "na")

fault_line(lon) = 0.3995299647473563 * lon + 46.31169212690954


grid = collect(Iterators.product(-82.:0.2:-79., 12.:0.2:16.))
glons = vec([g[1] for g in grid])
glats = vec([g[2] for g in grid])

partials = Oiler.Elastic.calc_locking_effects_per_fault(ss, glons, glats);

pv = [pole.x; pole.y; pole.z];

pred_vecs = [part * pv for part in partials];

pe, pn = [v[1] for v in pred_vecs], [v[2] for v in pred_vecs]

ca_inds = [i for (i, glon) in enumerate(glons) if fault_line(glon) > glats[i] ]

ca_lons = glons[ca_inds]
ca_lats = glats[ca_inds]

ca_vels = Oiler.predict_block_vels(ca_lons, ca_lats, pole)
ca_ve = [v.ve for v in ca_vels]
ca_vn = [v.vn for v in ca_vels]

ve, vn = pe .* 1., pn .* 1.

for (i, ci) in enumerate(ca_inds)
    ve[ci] += ca_ve[i]
    vn[ci] += ca_vn[i]
end

# results in meters per year

figure()
#quiver(glons, glats, pe, pn)
quiver(glons, glats, ve, vn)
plot(strike_slip_trace[:,1], strike_slip_trace[:,2])
axis("equal")
show()
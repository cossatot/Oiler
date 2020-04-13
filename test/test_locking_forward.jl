using Revise

using Oiler

using PyPlot


pole = Oiler.PoleCart(x = -4.000480063406787e-10, y = -3.686194245405466e-9,
                      z = 3.100995264274769e-9, fix = "na", mov = "ca");

thrust_trace = [-81.905471 12.945880; -81.554797 12.114126];
thrust = Oiler.Fault(trace = thrust_trace, dip_dir = "SW", dip = 20., hw = "na",
                      fw = "ca")

v01 = Oiler.VelocityVectorSphere(lond = -83.008, latd = 12.019, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "01", vel_type = "gnss");
v02 = Oiler.VelocityVectorSphere(lond = -82.157, latd = 12.360, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "02", vel_type = "gnss");
v03 = Oiler.VelocityVectorSphere(lond = -81.901, latd = 12.462, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "03", vel_type = "gnss");
v04 = Oiler.VelocityVectorSphere(lond = -81.815, latd = 12.496, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "04", vel_type = "gnss");
v05 = Oiler.VelocityVectorSphere(lond = -81.739, latd = 12.527, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "05", vel_type = "gnss");
v06 = Oiler.VelocityVectorSphere(lond = -81.721, latd = 12.533, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "06", vel_type = "gnss");
v07 = Oiler.VelocityVectorSphere(lond = -81.645, latd = 12.564, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "07", vel_type = "gnss");
v08 = Oiler.VelocityVectorSphere(lond = -81.559, latd = 12.598, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "08", vel_type = "gnss");
v09 = Oiler.VelocityVectorSphere(lond = -81.304, latd = 12.701, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "09", vel_type = "gnss");
v10 = Oiler.VelocityVectorSphere(lond = -80.450, latd = 13.041, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "10", vel_type = "gnss");

gnss_vels = [v01, v02, v03, v04, v05, v06, v07, v08, v09, v10];

glons, glats = Oiler.get_coords_from_vel_array(gnss_vels)

partials = Oiler.Elastic.calc_locking_effects_per_fault(thrust, glons, glats);

pv = [pole.x; pole.y; pole.z];

pred_vecs = [part * pv for part in partials];

pe, pn = [v[1] for v in pred_vecs], [v[2] for v in pred_vecs]

figure()
quiver(glons, glats, pe, pn, scale = 100)
plot(thrust_trace[:,1], thrust_trace[:,2])
axis("equal")
show()
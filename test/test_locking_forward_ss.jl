using Revise

using Oiler

using PyPlot


pole = Oiler.PoleCart(x = -4.000480063406787e-10, y = -3.686194245405466e-9,
                      z = 3.100995264274769e-9, fix = "na", mov = "ca");

# thrust_trace = [-81.905471 12.945880; -81.554797 12.114126];
# thrust = Oiler.Fault(trace = thrust_trace, dip_dir = "SW", dip = 20., hw = "na",
#                      fw = "ca")

strike_slip_trace = [-81.5 11.75; -80.650 12.091; -79.798 12.430];
ss = Oiler.Fault(trace = strike_slip_trace, dip_dir = "N", dip = 89., hw = "ca", fw = "na")

v11 = Oiler.VelocityVectorSphere(lond = -81.178, latd = 13.338, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "11", vel_type = "gnss");
v12 = Oiler.VelocityVectorSphere(lond = -80.825, latd = 12.507, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "12", vel_type = "gnss");
v13 = Oiler.VelocityVectorSphere(lond = -80.720, latd = 12.257, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "13", vel_type = "gnss");
v14 = Oiler.VelocityVectorSphere(lond = -80.685, latd = 12.174, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "14", vel_type = "gnss");
v15 = Oiler.VelocityVectorSphere(lond = -80.653, latd = 12.099, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "15", vel_type = "gnss");
v16 = Oiler.VelocityVectorSphere(lond = -80.646, latd = 12.083, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "16", vel_type = "gnss");
v17 = Oiler.VelocityVectorSphere(lond = -80.615, latd = 12.008, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "17", vel_type = "gnss");
v18 = Oiler.VelocityVectorSphere(lond = -80.580, latd = 11.925, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "18", vel_type = "gnss");
v19 = Oiler.VelocityVectorSphere(lond = -80.475, latd = 11.675, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "19", vel_type = "gnss");
v20 = Oiler.VelocityVectorSphere(lond = -80.126, latd = 10.843, ve = 0., 
    vn = 0., fix = "na", mov = "ca", name = "20", vel_type = "gnss");

gnss_vels = [v11, v12, v13, v14, v15, v16, v17, v18, v19, v20];

PvGbs = [Oiler.BlockRotations.build_PvGb_vel(v) for v in gnss_vels];

PvGbf = Oiler.BlockRotations.build_PvGb_vel(
    Oiler.fault_to_vel(ss)
);


Pf = Oiler.Faults.build_velocity_projection_matrix(ss.strike, ss.dip)
# Pf = [1. 0. 0; 0. 1. 0.; 0. 0. 1.];
Pa = Oiler.Faults.build_strike_rot_matrix(ss.strike);

glons, glats = Oiler.get_coords_from_vel_array(gnss_vels)

partials = Oiler.Elastic.calc_locking_effects_per_fault(ss, glons, glats);

pv = [pole.x; pole.y; pole.z];

pred_vecs = [part * pv for part in partials];

pe, pn = [v[1] for v in pred_vecs], [v[2] for v in pred_vecs]

figure()
quiver(glons, glats, pe, pn, scale = 100)
plot(strike_slip_trace[:,1], strike_slip_trace[:,2])
axis("equal")
show()
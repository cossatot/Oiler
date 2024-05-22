using Test

using PyPlot

using Oiler

"""
This script uses several Euler poles 
that form good velocity triangles, to test various aspects of the
pole inversion capabilities of Oiler.

The sites are approximately located on the plate boundaries of interest,
assuming no division of India into India+Capricorn, and Africa into
Nubia+Somalia+whatever.

The tests essentially take the poles, predict velocities at each site
based on them, and then try to invert this data to recover the poles.
"""


af_an_s = Oiler.PoleSphere(lat = 3.661, lon = -31.669, rotrate = 0.158, mov = "an", 
                                fix = "af");
in_af_s = Oiler.PoleSphere(lat = 18.915, lon = 47.060, rotrate = 0.606, mov = "af", 
                                fix = "in");
in_an_s = Oiler.PoleSphere(lat = 18.328, lon = 32.738, rotrate = 0.657, mov = "an",
                                fix = "in");
in_ar_s = Oiler.PoleSphere(lat = -7.484, lon = -174.932, rotrate = 0.356, mov = "ar",
                                fix = "in");
ar_af_s = Oiler.PoleSphere(lon = 31.290, lat = 15.625, rotrate = 0.901, mov = "af", 
                                fix = "ar");
eu_ar_s = Oiler.PoleSphere(lon = 64.748, lat = 26.733, rotrate = 0.368, mov = "ar",
                                fix = "eu");
eu_in_s = Oiler.PoleSphere(lon = 33.775, lat = 19.727, rotrate = 0.628, mov = "in", 
                                fix = "eu");
# af_eu_s = Oiler.PoleSphere(lat = 25.3, lon = -21.2, rotrate = 0.10, mov = "af", 
#                                fix = "eu");


af_an_sites = [13.98 -52.17;
               18.15 -52.45;
               25.56 -52.82;
               30.58 -49.12;
               52.94 -36.14;
               64.49 -27.83];

af_an_vels = Oiler.predict_block_vels(af_an_sites[:,1], af_an_sites[:,2], af_an_s);


in_af_sites = [69.61 -24.31;
               68.90 -19.79;
               65.51 -18.85;
               67.33 -16.72;
               67.70 -8.79;
               66.97 -0.17;
               62.15 4.79;
               57.84 11.81];

               
in_af_vels = Oiler.predict_block_vels(in_af_sites[:,1], in_af_sites[:,2],
                                      in_af_s);


in_an_sites = [130.01 -50.55;
               104.13 -48.55;
               84.72 -41.93;
               75.65 -30.48;
               70.14 -25.81];

in_an_vels = Oiler.predict_block_vels(in_an_sites[:,1], in_an_sites[:,2],
                                      in_an_s);


in_ar_sites = [61.05 16.28;
               58.85 14.90;
               59.99 16.33;
               60.93 19.72;
               62.00 21.30;
               64.36 23.77];

in_ar_vels = Oiler.predict_block_vels(in_ar_sites[:,1], in_ar_sites[:,2],
                                      in_ar_s);


ar_af_sites = [54.74 14.81;
               51.52 13.58;
               46.17 12.10;
               44.46 11.98;
               42.21 14.81;
               40.03 18.08
               35.52 29.20];

ar_af_vels = Oiler.predict_block_vels(ar_af_sites[:,1], ar_af_sites[:,2],
                                      ar_af_s);
eu_ar_sites = [38.31 38.05;
               41.70 38.18;
               46.04 36.68;
               53.50 27.01;
               61.47 24.11];

eu_ar_vels = Oiler.predict_block_vels(eu_ar_sites[:,1], eu_ar_sites[:,2],
                                      eu_ar_s);

eu_in_sites = [79.33 29.0; 
         82.55 27.67;
         75.89 32.08;
         94.79 27.85;
         72.78 35.27;
         79.1  46.03];

eu_in_vels = Oiler.predict_block_vels(eu_in_sites[:,1], eu_in_sites[:,2],
                                      eu_in_s);


vels = reduce(vcat, (af_an_vels, 
                     in_af_vels, 
                     in_an_vels, 
                     in_ar_vels,
                     ar_af_vels,
                     eu_ar_vels,
                     eu_in_vels));

vel_groups = Oiler.group_vels_by_fix_mov(vels);

gg = Oiler.make_block_PvGb_from_vels(vel_groups);

poles = Oiler.solve_block_invs_from_vel_groups(vel_groups);

af_an_pred = Oiler.pole_cart_to_sphere(poles[("af", "an")]);
in_af_pred = Oiler.pole_cart_to_sphere(poles[("in", "af")]);
in_an_pred = Oiler.pole_cart_to_sphere(poles[("in", "an")]);


ps = [Oiler.pole_cart_to_sphere(p) for (k, p) in poles];


ve_pred, vn_pred = Oiler.predict_vels_from_poles(gg, [p for (k, p) in poles]);

ve_obs, vn_obs = Oiler.predict_vels_from_poles(gg, [af_an_s, in_af_s, in_an_s, 
                                                    in_ar_s, ar_af_s, eu_ar_s, 
                                                    eu_in_s]);

vel_lons = reduce(vcat, [[v.lon for v in vel_groups[key]] for key in gg["keys"]])
vel_lats = reduce(vcat, [[v.lat for v in vel_groups[key]] for key in
gg["keys"]])

figure()
quiver(vel_lons, vel_lats, ve_pred, vn_pred, scale = 300., color = "red")
quiver(vel_lons, vel_lats, ve_obs, vn_obs, scale = 300., color = "blue")
show()


function compare_poles(p1, p2)
    @test round(p1.lon; digits = 1) == round(p2.lon; digits = 1)
    @test round(p1.lat; digits = 1) == round(p2.lat; digits = 1)
    @test round(p1.rotrate; digits = 2) == round(p2.rotrate; digits = 2)
end

@testset "five_plates.jl comparing input vs. output poles" begin
compare_poles(af_an_s, af_an_pred)
compare_poles(in_af_s, in_af_pred)
compare_poles(in_an_s, in_an_pred)
end
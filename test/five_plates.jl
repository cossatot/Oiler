using Oiler

"""
This script uses several Euler poles from Cox and Hart (1986),
intially derived from Minster and Jordan (1978),
that form a good velocity triangle, to test various aspects of the
pole inversion capabilities of Oiler.

The sites are approximately located on the plate boundaries of interest,
assuming no division of India into India+Capricorn, and Africa into
Nubia+Somalia+whatever.

The tests essentially take the poles, predict velocities at each site
based on them, and then try to invert this data to recover the poles.
"""


an_af_s = Oiler.EulerPoleSphere(latd=9.6, lond=-41.6, rotrate=0.15, mov="an", 
                                fix="af");
af_in_s = Oiler.EulerPoleSphere(latd=17.3, lond=46., rotrate=0.64, mov="af", 
                                fix="in");
an_in_s = Oiler.EulerPoleSphere(latd=18.7, lond=32.5, rotrate=0.67, mov="an",
                                fix="in");
ar_in_s = Oiler.EulerPoleSphere(latd=-16.3, lond=-152.2, rotrate=0.61, mov="ar",
                                fix="in");
ar_af_s = Oiler.EulerPoleSphere(lond=-61.5, latd=-5.6, rotrate=0.19, mov="af", 
                                fix="ar");
af_eu_s = Oiler.EulerPoleSphere(latd=25.3, lond=-21.2, rotrate=0.10, mov="af", 
                                fix="eu");
ar_eu_s = Oiler.EulerPoleSphere(lond= 91.79, latd=25.17, rotrate=0.15, mov="ar",
                                fix="eu");
eu_in_s = Oiler.EulerPoleSphere(lond=38.4, latd=19.7, rotrate=0.70, mov="in", 
                                fix="eu");


an_af_sites = [13.98 -52.17;
               18.15 -52.45;
               25.56 -52.82;
               30.58 -49.12;
               52.94 -36.14;
               64.49 -27.83];

an_af_vels = Oiler.predict_block_vels(an_af_sites[:,1], an_af_sites[:,2], an_af_s);


af_in_sites = [69.61 -24.31;
               68.90 -19.79;
               65.51 -18.85;
               67.33 -16.72;
               67.70 -8.79;
               66.97 -0.17;
               62.15 4.79;
               57.84 11.81];
               
af_in_vels = Oiler.predict_block_vels(af_in_sites[:,1], af_in_sites[:,2],
                                      af_in_s);


an_in_sites = [130.01 -50.55;
               104.13 -48.55;
               84.72 -41.93;
               75.65 -30.48;
               70.14 -25.81];

an_in_vels = Oiler.predict_block_vels(an_in_sites[:,1], an_in_sites[:,2],
                                      an_in_s);


ar_in_sites = [61.05 16.28;
               58.85 14.90;
               59.99 16.33;
               60.93 19.72;
               62.00 21.30;
               64.36 23.77];

ar_in_vels = Oiler.predict_block_vels(ar_in_sites[:,1], ar_in_sites[:,2],
                                      ar_in_s);


ar_af_sites = [54.74 14.81;
               51.52 13.58;
               46.17 12.10;
               44.46 11.98;
               42.21 14.81;
               40.03 18.08
               35.52 29.20];

ar_af_vels = Oiler.predict_block_vels(ar_af_sites[:,1], ar_af_sites[:,2],
                                      ar_af_s);
ar_eu_sites = [38.31 38.05;
               41.70 38.18;
               46.04 36.68;
               53.50 27.01;
               61.47 24.11];

ar_eu_vels = Oiler.predict_block_vels(ar_eu_sites[:,1], ar_eu_sites[:,2],
                                      ar_eu_s);

eu_in_sites = [79.33 29.0; 
         82.55 27.67;
         75.89 32.08;
         94.79 27.85;
         72.78 35.27;
         79.1  46.03];

eu_in_vels = Oiler.predict_block_vels(eu_in_sites[:,1], eu_in_sites[:,2],
                                      eu_in_s);


vels = reduce(vcat, (an_af_vels, af_in_vels, an_in_vels));#, ar_in_vels));#, ar_af_vels,
                     #ar_eu_vels, eu_in_vels));

vel_groups = Oiler.group_vels_by_fix_mov(vels);

poles = Oiler.solve_block_invs_from_vel_groups(vel_groups);

an_af_pred = Oiler.euler_pole_cart_to_sphere(poles[("af", "an")]);
af_in_pred = Oiler.euler_pole_cart_to_sphere(poles[("in", "af")]);
an_in_pred = Oiler.euler_pole_cart_to_sphere(poles[("in", "an")]);


function compare_poles(p1, p2)
    if round(p1.lond; digits=1) != round(p2.lond; digits=1)
        print("aaah lons", p1.fix, p1.mov)
    end
    if round(p1.latd; digits=1) != round(p2.latd; digits=1)
        print("aaah lats", p1.fix, p1.mov)
    end
    if round(p1.rotrate; digits=2) != round(p2.rotrate; digits=2)
        print("aaah rate", p1.fix, p1.mov)
    end
end

compare_poles(an_af_s, an_af_pred);
compare_poles(af_in_s, af_in_pred);
compare_poles(an_in_s, an_in_pred);

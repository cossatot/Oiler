using Oiler
using PyPlot


ind_eur_pole = Oiler.PoleSphere(lon=13.9, lat=26.5, rotrate=0.341, fix="eur", mov="ind");


sites = [79.33 29.0; 
         82.55 27.67;
         75.89 32.08;
         94.79 27.85;
         72.78 35.27;
         79.1  46.03];

site_lons = sites[:,1];
site_lats = sites[:,2];

pred_vels = Oiler.predict_block_vels(site_lons, site_lats, ind_eur_pole);

ve_pred = [pv.ve for pv in pred_vels];
vn_pred = [pv.vn for pv in pred_vels];

PvGb = Oiler.build_PvGb_from_vels(pred_vels);

V = Oiler.build_vel_column_from_vels(pred_vels);

oh = PvGb \ V;

pred_pole = Oiler.pole_cart_to_sphere(Oiler.PoleCart(x=oh[1], y=oh[2], z=oh[3]));


figure()
quiver(site_lons, site_lats, ve_pred, vn_pred)
show()
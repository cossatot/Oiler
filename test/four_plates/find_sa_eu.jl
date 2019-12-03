using Revise

using DataFrames
using CSV
using PyPlot

using Oiler

vels = CSV.read("./data/fault_vels.csv");

function vel_from_row(row::DataFrameRow)
    Oiler.VelocityVectorSph(lond = row.lon, latd = row.lat, 
                            ve = row.ve, vn = row.vn)
end

af_eu_vels = Array{Oiler.VelocityVectorSph,1}();
na_af_vels = Array{Oiler.VelocityVectorSph,1}();
na_eu_vels = Array{Oiler.VelocityVectorSph,1}();
na_sa_vels = Array{Oiler.VelocityVectorSph,1}();
af_sa_vels = Array{Oiler.VelocityVectorSph,1}();

for i in 1:size(vels, 1)
    row = vels[i,:]
    if row.plates == "af_eu"
        push!(af_eu_vels, vel_from_row(row));
    elseif row.plates == "na_af"
        push!(na_af_vels, vel_from_row(row));
    elseif row.plates == "na_eu"
        push!(na_eu_vels, vel_from_row(row));
    elseif row.plates == "na_sa"
        push!(na_sa_vels, vel_from_row(row));
    elseif row.plates == "af_sa"
        push!(af_sa_vels, vel_from_row(row));
    end
end

af_eu_PvGb = Oiler.build_PvGb_from_vels(af_eu_vels);
na_af_PvGb = Oiler.build_PvGb_from_vels(na_af_vels);
na_eu_PvGb = Oiler.build_PvGb_from_vels(na_eu_vels);
na_sa_PvGb = Oiler.build_PvGb_from_vels(na_sa_vels);
af_sa_PvGb = Oiler.build_PvGb_from_vels(af_sa_vels);



big_PvGb = Oiler.diagonalize_matrices((af_eu_PvGb, na_af_PvGb, na_eu_PvGb, na_sa_PvGb,
                                 af_sa_PvGb));


big_vels = Oiler.build_vel_column_from_vels([af_eu_vels; na_af_vels; na_eu_vels;
                                             na_sa_vels; af_sa_vels]);

omegas = big_PvGb \ big_vels;

af_eu_pole = Oiler.EulerPoleCart(x = omegas[1], y = omegas[2], z = omegas[3],
                                fix = "af", rel = "eu");
na_af_pole = Oiler.EulerPoleCart(x = omegas[4], y = omegas[5], z = omegas[6],
                                 fix = "na", rel = "af");
na_eu_pole = Oiler.EulerPoleCart(x = omegas[7], y = omegas[8], z = omegas[9],
                                 fix = "na", rel = "eu");
na_sa_pole = Oiler.EulerPoleCart(x = omegas[10], y = omegas[11], z = omegas[12],
                                 fix = "na", rel = "sa");
af_sa_pole = Oiler.EulerPoleCart(x = omegas[13], y = omegas[14], z = omegas[15],
                                 fix = "af", rel = "sa");

sa_eu_pole_1 = Oiler.euler_pole_cart_to_sphere(Oiler.add_poles(na_sa_pole, na_eu_pole));
sa_eu_pole_2 = Oiler.euler_pole_cart_to_sphere(Oiler.add_poles(na_af_pole, af_eu_pole));

na_eu_lats = [vel.latd for vel in na_eu_vels];
na_eu_lons = [vel.lond for vel in na_eu_vels];

na_eu_pred = Oiler.predict_block_vels(na_eu_lons, na_eu_lats, na_eu_pole);
na_eu_pole_af = Oiler.add_poles(na_af_pole, af_eu_pole)
na_eu_pred_af = Oiler.predict_block_vels(na_eu_lons, na_eu_lats, na_eu_pole_af);


figure()
quiver(na_eu_lons, na_eu_lats, [v.ve for v in na_eu_vels], 
       [v.vn for v in na_eu_vels], color = "blue", scale = 100)
quiver(na_eu_lons, na_eu_lats, [v.ve for v in na_eu_pred], 
       [v.vn for v in na_eu_pred], color = "red", scale = 100)
quiver(na_eu_lons, na_eu_lats, [v.ve for v in na_eu_pred_af], 
       [v.vn for v in na_eu_pred_af], color = "green", scale = 100)
show()
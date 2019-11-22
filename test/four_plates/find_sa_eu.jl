using DataFrames
using CSV
using PyPlot

using Oiler

vels = CSV.read("./data/fault_vels.csv");

function vel_from_row(row::DataFrameRow)
    Oiler.VelocityVectorSph(row.lon, row.lat, row.ve, row.vn, 0., 0., 0., 0., )
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

function diagonalize_matrices(matrices)

    rowz = [size(m, 1) for m in matrices]
    colz = [size(m, 2) for m in matrices]

    n_rows = sum(rowz)
    n_cols = sum(colz)

    big_mat = zeros(n_rows, n_cols)

    i_row = 1
    i_col = 1

    for (i, mat) in enumerate(matrices)
        row_end = i_row + rowz[i] - 1
        col_end = i_col + colz[i] - 1

        big_mat[i_row:row_end, i_col:col_end] = mat

        i_row = row_end + 1
        i_col = col_end + 1
    end
    return big_mat
end


big_PvGb = diagonalize_matrices((af_eu_PvGb, na_af_PvGb, na_eu_PvGb, na_sa_PvGb,
                                 af_sa_PvGb));


big_vels = Oiler.build_vel_column_from_vels([af_eu_vels; na_af_vels; na_eu_vels;
                                             na_sa_vels; af_sa_vels]);

omegas = big_PvGb \ big_vels;

af_eu_pole = Oiler.EulerPoleCart(omegas[1], omegas[2], omegas[3]);
na_af_pole = Oiler.EulerPoleCart(omegas[4], omegas[5], omegas[6]);
na_eu_pole = Oiler.EulerPoleCart(omegas[7], omegas[8], omegas[9]);
na_sa_pole = Oiler.EulerPoleCart(omegas[10], omegas[11], omegas[12]);
af_sa_pole = Oiler.EulerPoleCart(omegas[13], omegas[14], omegas[15]);

sa_eu_pole_1 = Oiler.euler_pole_cart_to_sphere(Oiler.add_poles(na_sa_pole, na_eu_pole));
sa_eu_pole_2 = Oiler.euler_pole_cart_to_sphere(Oiler.add_poles(na_af_pole, af_eu_pole));

sa_eu_vels = CSV.read(("./data/sa_eu_vels.csv"));
sa_eu_gps = [vel_from_row(sa_eu_vels[i,:]) for i in 1:size(sa_eu_vels, 1)];

sa_eu_v_lons = [v.lond for v in sa_eu_gps];
sa_eu_v_lats = [v.latd for v in sa_eu_gps];

sa_eu_1_pred = Oiler.predict_block_vels(sa_eu_v_lons, sa_eu_v_lats, 
                                        sa_eu_pole_1);

sa_eu_2_pred = Oiler.predict_block_vels(sa_eu_v_lons, sa_eu_v_lats, 
                                        sa_eu_pole_2);

sa_eu_1p_ve = [pv.ve for pv in sa_eu_1_pred];
sa_eu_1p_vn = [pv.vn for pv in sa_eu_1_pred];

sa_eu_2p_ve = [pv.ve for pv in sa_eu_2_pred];
sa_eu_2p_vn = [pv.vn for pv in sa_eu_2_pred];

v_scale = 200.

figure()
quiver(sa_eu_v_lons, sa_eu_v_lats, [v.ve for v in sa_eu_gps], 
       [v.vn for v in sa_eu_gps], color = "black", scale = v_scale)
quiver(sa_eu_v_lons, sa_eu_v_lats, sa_eu_1p_ve, sa_eu_1p_vn, color = "red", scale = v_scale)
quiver(sa_eu_v_lons, sa_eu_v_lats, sa_eu_2p_ve, sa_eu_2p_vn, color = "blue",
scale = v_scale)
show()
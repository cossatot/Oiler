include("..\\src\\BlockRotations.jl")
#using .BlockRotations
using PyPlot

struct Vec
    latd::AbstractFloat
    lond::AbstractFloat
    ve::AbstractFloat
    vn::AbstractFloat
    vu::AbstractFloat
    ee::AbstractFloat
    en::AbstractFloat
    eu::AbstractFloat
end

AVES = Vec(15.67, -63.62, 13.2, 14.2, 0., 0., 0., 0.);
ROJO = Vec(17.9, -71.67, 13.5, 6.9, 0., 0., 0., 0.);
SANA = Vec(12.53, -81.73, 14.8, 7.4, 0., 0., 0., 0.);
CRO1 = Vec(17.76, -64.58, 9.8, 13.3, 0., 0., 0., 0.);

pv_aves = build_Pv_deg(AVES.lond, AVES.latd);
pv_rojo = build_Pv_deg(ROJO.lond, ROJO.latd);
pv_sana = build_Pv_deg(SANA.lond, SANA.latd);
pv_cro1 = build_Pv_deg(CRO1.lond, CRO1.latd);

gb_aves = build_Gb_deg(AVES.lond, AVES.latd);
gb_rojo = build_Gb_deg(ROJO.lond, ROJO.latd);
gb_sana = build_Gb_deg(SANA.lond, SANA.latd);
gb_cro1 = build_Gb_deg(CRO1.lond, CRO1.latd);


A = [pv_aves * gb_aves;
     pv_rojo * gb_rojo;
     pv_sana * gb_sana;
     pv_cro1 * gb_cro1];

V = [AVES.ve; AVES.vn; AVES.vu;
     ROJO.ve; ROJO.vn; ROJO.vu;
     SANA.ve; SANA.vn; SANA.vu
     CRO1.ve; CRO1.vn; CRO1.vu];

#V *= 1e3;

omega_hat = A \ V;

rrate = sqrt(sum(omega_hat.^2));
omega_hat_unit = omega_hat ./ rrate;
#om_lon = atand(omega_hat_unit[2],omega_hat_unit[1]);
#om_lat = atand(omega_hat_unit[3] ./ sqrt(omega_hat_unit[1]^2 + omega_hat_unit[2]^2));

om_lon = atand(omega_hat_unit[1], omega_hat_unit[2]) + 360.;
om_lat = acosd(omega_hat_unit[3]);

println((om_lon, om_lat, rrate))

xs = [AVES.lond ROJO.lond SANA.lond CRO1.lond];
ys = [AVES.latd ROJO.latd SANA.latd CRO1.latd];

ve_obs = [AVES.ve ROJO.ve SANA.ve CRO1.ve];
vn_obs = [AVES.vn ROJO.vn SANA.vn CRO1.vn];

ve_pred = [(pv_aves * gb_aves * omega_hat)[1]
           (pv_rojo * gb_rojo * omega_hat)[1]
           (pv_sana * gb_sana * omega_hat)[1]
           (pv_cro1 * gb_cro1 * omega_hat)[1]];

vn_pred = [(pv_aves * gb_aves * omega_hat)[2]
           (pv_rojo * gb_rojo * omega_hat)[2]
           (pv_sana * gb_sana * omega_hat)[2]
           (pv_cro1 * gb_cro1 * omega_hat)[2]];


figure()
quiver(xs, ys, ve_pred, vn_pred, color="blue", scale=100.)
quiver(xs, ys, ve_obs, vn_obs, color="red", scale=100.)
show()
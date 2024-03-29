using PyPlot

using CSV
using JSON
using DataFrames

using Oiler


pole = Oiler.PoleCart(x=-4.000480063406787e-10, y=-3.686194245405466e-9,
                      z=3.100995264274769e-9, fix="na", mov="ca");

pv = [pole.x; pole.y; pole.z];

# load tris
tri_json = JSON.parsefile("./test_data/fake_na_ca/ca_na_fake_la_tris.geojson",
    dicttype=Dict)

function tri_from_geojson(tri_feature)
    coords = tri_feature["geometry"]["coordinates"][1]
    Oiler.Tris.Tri(p1=Float64.(coords[1]), 
                   p2=Float64.(coords[2]), 
                   p3=Float64.(coords[3]), 
                   name=tri_feature["properties"]["fid"])
end

tris = map(tri_from_geojson, tri_json["features"])

# load vel points
vels = CSV.read("./test_data/fake_na_ca/ca_na_fake_pts.csv", DataFrame)
vlon = [vels[i,:].X for i in 1:size(vels, 1)]
vlat = [vels[i,:].Y for i in 1:size(vels, 1)]

ca_idx = [i for i in 1:size(vels, 1) if vels[i,4] == "ca"]

# calc plate vel at each point
pe = zeros(size(vlon))
pn = zeros(size(vlon))

ca_vels = Oiler.predict_block_vels(vlon[ca_idx], vlat[ca_idx], pole)
ca_ve = [v.ve for v in ca_vels]
ca_vn = [v.vn for v in ca_vels]

for (i, ci) in enumerate(ca_idx)
    pe[ci] += ca_ve[i]
    pn[ci] += ca_vn[i]
end

fault_err_scale = 1.0e0

# load faults
ss1 = Oiler.Fault(trace=[-81.5 17.76; -80.650 18.091; -79.798 18.430],
    dip_dir="S", dip=89., hw="ca", fw="na", name="ss1",
    dextral_rate=-20., dextral_err=20. * fault_err_scale, 
    extension_rate=0., extension_err=5. * fault_err_scale)
ss2 = Oiler.Fault(trace=[-83.5 17.76; -82.5  17.75; -81.5 17.75],
    dip_dir="S", dip=89., hw="ca", fw="na", name="ss2",
    dextral_rate=-20., dextral_err=20. * fault_err_scale, 
    extension_rate=0., extension_err=5. * fault_err_scale)
# th1 = Oiler.Fault(trace=[ -79.798 18.430; -78.08 14.72],
th1 = Oiler.Fault(trace=[ -79.798 18.430; -78.76 17.54; -78.08 14.72],
    dip_dir="W", dip=20., hw="ca", fw="na", name="th1",
    dextral_rate=0., dextral_err=10. * fault_err_scale, 
    extension_rate=-20., extension_err=20. * fault_err_scale)

# redo with purrfect rates
ss1_rl, ss1_ex = Oiler.Faults.get_fault_slip_rate_from_pole(ss1, -pole)
ss2_rl, ss2_ex = Oiler.Faults.get_fault_slip_rate_from_pole(ss2, pole)
th1_rl, th1_ex = Oiler.Faults.get_fault_slip_rate_from_pole(th1, pole)

ss1 = Oiler.Fault(trace=[-81.5 17.76; -80.650 18.091; -79.798 18.430],
    dip_dir="S", dip=89., hw="ca", fw="na", name="ss1",
    dextral_rate=ss1_rl, dextral_err=20. * fault_err_scale, 
    extension_rate=ss1_ex, extension_err=5. * fault_err_scale)
ss2 = Oiler.Fault(trace=[-83.5 17.76; -82.5  17.75; -81.5 17.75],
    dip_dir="S", dip=89., hw="ca", fw="na", name="ss2",
    dextral_rate=ss2_rl, dextral_err=20. * fault_err_scale, 
    extension_rate=ss2_ex, extension_err=5. * fault_err_scale)
# th1 = Oiler.Fault(trace=[ -79.798 18.430; -78.08 14.72],
th1 = Oiler.Fault(trace=[ -79.798 18.430; -78.76 17.54; -78.08 14.72],
    dip_dir="W", dip=20., hw="ca", fw="na", name="th1",
    dextral_rate=th1_rl, dextral_err=10. * fault_err_scale, 
    extension_rate=th1_ex, extension_err=20. * fault_err_scale)


faults = [ss1; ss2; th1]
# calc partials for each fault

ss1_part = Oiler.Elastic.calc_locking_effects_per_fault(ss1, vlon, vlat)
ss2_part = Oiler.Elastic.calc_locking_effects_per_fault(ss2, vlon, vlat)
th1_part = Oiler.Elastic.calc_locking_effects_segmented_fault(th1, vlon, vlat)

# add together
parts = ss1_part + ss2_part + th1_part
lock_vels = [part * pv for part in parts]
le, ln = [v[1] for v in lock_vels], [v[2] for v in lock_vels]

ve = pe - le
vn = pn - ln


# inversion
final_vels = [Oiler.VelocityVectorSphere(lon=lon, lat=vlat[i],  
                                         ve=ve[i], vn=vn[i],
                                         # ve=pe[i], vn=pn[i],
                                         ee=0.01, en=0.01, # eu=0.001,
                                         fix="fix", 
                                         mov=(if i in ca_idx; "ca" else "na" end),
                                         # mov="ca", 
                                         vel_type="GNSS", name=string(i))
               for (i, lon) in enumerate(vlon)]

fault_vels = [Oiler.fault_to_vel(fault) for fault in (ss1, ss2, th1)]

all_vels = vcat(final_vels, fault_vels)

vel_groups = Oiler.group_vels_by_fix_mov(all_vels)

results = Oiler.solve_block_invs_from_vel_groups(vel_groups, 
    # faults=faults, 
    faults=[ss1, ss2],
    tris=tris,
    tri_distance_weight=0.1,
    regularize_tris=false,
    weighted=true,
    predict_vels=true,
    pred_se=true)

in_g_vels = Oiler.Utils.make_df_from_vel_array(final_vels)
pred_g_vels = Oiler.Utils.make_gnss_df_from_vel_groups(results["predicted_vels"])



figure()
quiver(vlon, vlat, ve, vn, scale=300, color="black")
# quiver(vlon, vlat, pve, pvn, scale=300, color="red", alpha=0.5)
quiver(pred_g_vels.lon, pred_g_vels.lat, pred_g_vels.ve, pred_g_vels.vn,
       scale=300, color="red", alpha=0.5)
plot(ss1.trace[:,1], ss1.trace[:,2])
plot(ss2.trace[:,1], ss2.trace[:,2])
plot(th1.trace[:,1], th1.trace[:,2])
axis("equal")
show()
display(gcf())

Oiler.pole_cart_to_sphere(results["poles"][("ca", "na")])

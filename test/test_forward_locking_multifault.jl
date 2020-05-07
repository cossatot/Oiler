using Revise
using PyPlot

using CSV
using DataFrames

using Oiler


pole = Oiler.PoleCart(x = -4.000480063406787e-10, y = -3.686194245405466e-9,
                      z = 3.100995264274769e-9, fix = "na", mov = "ca");

pv = [pole.x; pole.y; pole.z];


# load vel points
vels = CSV.read("./test_data/fake_na_ca/ca_na_fake_pts.csv")
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


# load faults
ss1 = Oiler.Fault(trace = [-81.5 17.76; -80.650 18.091; -79.798 18.430],
    dip_dir = "S", dip = 89., hw = "ca", fw = "na")
ss2 = Oiler.Fault(trace = [-83.5 17.76; -82.5  17.75; -81.5 17.75],
    dip_dir = "S", dip = 89., hw = "ca", fw = "na")
# th1 = Oiler.Fault(trace=[ -79.798 18.430; -78.08 14.72],
th1 = Oiler.Fault(trace = [ -79.798 18.430; -78.76 17.54; -78.08 14.72],
    dip_dir = "W", dip = 20., hw = "ca", fw = "na")

# calc partials for each fault

ss1_part = Oiler.Elastic.calc_locking_effects_per_fault(ss1, vlon, vlat)
ss2_part = Oiler.Elastic.calc_locking_effects_per_fault(ss2, vlon, vlat)
th1_part = Oiler.Elastic.calc_locking_effects_segmented_fault(th1, vlon, vlat)

# add together
parts = ss1_part + ss2_part + th1_part
lock_vels = [part * pv for part in parts]
le, ln = [v[1] for v in lock_vels], [v[2] for v in lock_vels]

ve = pe + le
vn = pn + ln

figure()
quiver(vlon, vlat, ve, vn)
plot(ss1.trace[:,1], ss1.trace[:,2])
plot(ss2.trace[:,1], ss2.trace[:,2])
plot(th1.trace[:,1], th1.trace[:,2])
axis("equal")
show()


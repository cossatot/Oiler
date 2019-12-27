using Revise

using DataFrames
using CSV
using PyPlot
using Random

using Oiler

Random.seed!(69)

const n_iters = 1000;

vels = Oiler.load_vels_from_csv("./data/fault_vels_err.csv");
vel_groups = Oiler.group_vels_by_fix_mov(vels);

poles = Oiler.solve_block_invs_from_vel_groups(vel_groups);

poles_iter = Oiler.solve_for_block_poles_iterative(vel_groups, n_iters);


gps_vels = Oiler.load_vels_from_csv("./data/eur_gps.csv");
gps_groups = Oiler.group_vels_by_fix_mov(gps_vels);

# make big pvgb
gp = Oiler.make_block_PvGb_from_vels(gps_groups)

pole_arr = [v for (k, v) in poles]

(vve, vvn) = Oiler.predict_vels_from_poles(gp, pole_arr)

ve_pred = zeros(length(gps_vels), n_iters)
vn_pred = zeros(length(gps_vels), n_iters)

for i in 1:n_iters
    (ve_pred[:,i], vn_pred[:,i]) = Oiler.predict_vels_from_poles(gp,
              poles_iter["poles"][i])
end

ve_obs = reduce(vcat, [[v.ve for v in gps_groups[key]] for key in gp["keys"]]);
vn_obs = reduce(vcat, [[v.vn for v in gps_groups[key]] for key in gp["keys"]]);

misfit(vep, veo, vnp, vno) = sum(sqrt.((veo .- vep).^2 .+ (vno .- vnp).^2))

scores = zeros(n_iters)

for i in 1:n_iters
    scores[i] = misfit(ve_pred[:,i], ve_obs, vn_pred[:,i], ve_obs)
end

best = argmin(scores)


gps_lons = reduce(vcat, [[v.lond for v in gps_groups[key]] for key in gp["keys"]])
gps_lats = reduce(vcat, [[v.latd for v in gps_groups[key]] for key in gp["keys"]])

 
figure()
quiver(gps_lons[1:10:end], gps_lats[1:10:end], ve_obs[1:10:end],
       vn_obs[1:10:end], scale = 100, color = "red")
quiver(gps_lons[1:10:end], gps_lats[1:10:end], ve_pred[1:10:end,best],
       vn_pred[1:10:end,best], scale = 100, color = "blue")
#quiver(gps_lons[1:10:end], gps_lats[1:10:end], vve[1:10:end],
#       vvn[1:10:end], scale = 100, color = "blue")
show()
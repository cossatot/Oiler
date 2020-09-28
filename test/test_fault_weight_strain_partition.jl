using Test

using PyPlot
# using Makie
using Setfield

using Oiler


# This script runs two near-identical plate motion inversions. The first



gnss_coords = [1. 1.;
               1. 1.5;
               1. 2.;
               1. 2.5;
               1. 3.;
               2. 1.;
               2. 1.5;
               2. 2.;
               2. 2.5;
               2. 3.;
               3. 1.;
               3. 2.;
               3. 3.];

pole = Oiler.PoleSphere(lon=-45., lat=45., rotrate=0.1, mov="fore", fix="ref")
pole_half = Oiler.PoleSphere(lon=-45., lat=45., rotrate=0.05, mov="sliv", fix="ref")


# foreland
pred_vels_fore = Oiler.predict_block_vels(gnss_coords[1:5,1], 
    gnss_coords[1:5,2], pole)

ppp = []
for vel in pred_vels_fore
    vel = @set vel.ee = 0.5
    vel = @set vel.en = 0.5
    push!(ppp, vel)
end

pred_vels_fore = ppp

# sliver
pred_vels_sliv = Oiler.predict_block_vels(gnss_coords[6:10,1], 
gnss_coords[6:10,2], pole_half)

ppp = []
for vel in pred_vels_sliv
    vel = @set vel.ee = 0.5
    vel = @set vel.en = 0.5
    push!(ppp, vel)
end

pred_vels_sliv = ppp

println(pred_vels_sliv[1])

# hinterland

v7 = Oiler.VelocityVectorSphere(
    lon=3., lat=1.,
    ve=0., vn=0.,
    ee=2., en=2.,
    fix="ref", mov="hint",
    vel_type="GNSS"
)

v8 = Oiler.VelocityVectorSphere(
    lon=3., lat=2.,
    ve=0., vn=0.,
    ee=2., en=2.,
    fix="ref", mov="hint",
    vel_type="GNSS"
)

v9 = Oiler.VelocityVectorSphere(
    lon=3., lat=3.,
    ve=0., vn=0.,
    ee=2., en=2.,
    fix="ref", mov="hint",
    vel_type="GNSS"
)


# faults

thrust_free = Oiler.Fault(
    trace=[1.5 0.01;
           1.5 4.],
    dip=10.,
    dip_dir="E",
    hw="sliv", fw="fore",
    name="thrust_free",
    dextral_rate=0.,
    dextral_err=50.,
    # extension_rate=0.,
    extension_rate=-7.67,
    extension_err=50.
)

thrust_fix = Oiler.Fault(
    trace=[1.5 0.01;
           1.5 4.],
    dip=10.,
    dip_dir="E",
    hw="sliv", fw="fore",
    name="thrust_fix",
    dextral_rate=0.,
    dextral_err=50.,
    extension_rate=-7.67,
    extension_err=0.
)

ss_free = Oiler.Fault(
    trace=[2.5 0.01;
           2.5 4.],
    dip=89.,
    dip_dir="W",
    hw="sliv", fw="hint",
    name="ss_free",
    dextral_rate=5.8,
    # dextral_rate=0.,
    dextral_err=50.,
    extension_rate=0.,
    extension_err=50.
)

ss_fix = Oiler.Fault(
    trace=[2.5 0.01;
           2.5 4.],
    dip=89.,
    dip_dir="W",
    hw="sliv", fw="hint",
    name="ss_fix",
    dextral_rate=5.8,
    dextral_err=0.,
    extension_rate=0.,
    extension_err=50.
)



function predict_gnss_vels(vel_groups, fault_groups, gnss_vels, poles)
    vg_keys = sort(collect(Tuple(keys(vel_groups))))

    vlon = [v.lon for v in gnss_vels]
    vlat = [v.lat for v in gnss_vels]

    locking_partial_groups = Dict()
    for vg in vg_keys
        if haskey(fault_groups, vg)
            locking_partial_groups[vg] = sum([
                Oiler.Elastic.calc_locking_effects_segmented_fault(fault, vlon, vlat)
                for fault in fault_groups[vg]
            ])
        # else
        #    println("Aint no fault here")
        end
    end
    
    pred_lock_vels_dict = Dict()
    
    for vg in keys(locking_partial_groups)
        pv = [poles[vg].x, poles[vg].y, poles[vg].z]
        pred_lock_vels_dict[vg] = [part * pv for part in
        locking_partial_groups[vg]]
    end
    
    pred_lock_vels = sum(values(pred_lock_vels_dict))
    
    
    ple, pln = [v[1] for v in pred_lock_vels], [v[2] for v in pred_lock_vels]
    
    gnss_pred_vels = []
    for vel in gnss_vels
        if haskey(poles, (vel.fix, vel.mov))
            push!(gnss_pred_vels, Oiler.Utils.get_vel_vec_from_pole(vel, poles[(vel.fix, vel.mov)]))
        else
            push!(gnss_pred_vels, Oiler.Utils.get_vel_vec_from_pole(vel, -poles[(vel.mov, vel.fix)]))
        end
    end
    
    pred_block_ve = [v[1] for v in gnss_pred_vels]
    pred_block_vn = [v[2] for v in gnss_pred_vels]
    
    pred_gnss_ve = ple + pred_block_ve
    pred_gnss_vn = ple + pred_block_vn

    (pred_gnss_ve, pred_gnss_vn)
end

# concat vels
gnss_vels = [pred_vels_fore[1], pred_vels_fore[2], pred_vels_fore[3],
             pred_vels_fore[4], pred_vels_fore[5],
             pred_vels_sliv[1], pred_vels_sliv[2], pred_vels_sliv[3],
             pred_vels_sliv[4], pred_vels_sliv[5],
             v7, v8, v9]
gnss_vels = convert(Array{VelocityVectorSphere}, gnss_vels)

fix_faults = [thrust_fix, ss_fix]
free_faults = [thrust_free, ss_free]

free_fault_vels = reduce(vcat, map(Oiler.fault_to_vels, free_faults))
fix_fault_vels = reduce(vcat, map(Oiler.fault_to_vels, fix_faults))

fix_vels = vcat(fix_fault_vels, gnss_vels)
free_vels = vcat(free_fault_vels, gnss_vels)

# group vels
fix_vel_groups = Oiler.group_vels_by_fix_mov(fix_vels)
free_vel_groups = Oiler.group_vels_by_fix_mov(free_vels)

vg_keys = sort(collect(Tuple(keys(fix_vel_groups))))

fault_groups_free = Oiler.Utils.group_faults(free_faults, vg_keys)
fault_groups_fix = Oiler.Utils.group_faults(fix_faults, vg_keys)


# solve
fix_poles = Oiler.solve_block_invs_from_vel_groups(fix_vel_groups,
    faults=fix_faults)

free_poles = Oiler.solve_block_invs_from_vel_groups(free_vel_groups; 
    faults=free_faults)

# find rates at faults
fix_rates = Oiler.Utils.get_fault_slip_rates_from_poles(fix_faults, fix_poles)
free_rates = Oiler.Utils.get_fault_slip_rates_from_poles(free_faults, free_poles)

# predict gnss vels w/ locking
gnss_pred_ve_fix, gnss_pred_vn_fix = predict_gnss_vels(fix_vel_groups, fault_groups_fix,
    gnss_vels, fix_poles)
gnss_pred_ve_free, gnss_pred_vn_free = predict_gnss_vels(free_vel_groups, fault_groups_free,
    gnss_vels, free_poles)


gnss_obs_ve = [v.ve for v in gnss_vels]
gnss_obs_vn = [v.vn for v in gnss_vels]

gnss_lon = [v.lon for v in gnss_vels]
gnss_lat = [v.lat for v in gnss_vels]


blocks_fix = Oiler.Solver.make_block_inversion_matrices_from_vels(fix_vel_groups)
blocks_free = Oiler.Solver.make_block_inversion_matrices_from_vels(free_vel_groups)


fix_rates


# @test ss_fix.dextral_rate 


figure(figsize=(10, 8))
 
plot(thrust_fix.trace[:,1], thrust_fix.trace[:,2])
plot(ss_fix.trace[:,1], ss_fix.trace[:,2])

quiver(gnss_lon, gnss_lat, gnss_obs_ve, gnss_obs_vn, color="k", scale=100)
quiver(gnss_lon, gnss_lat, gnss_pred_ve_free, gnss_pred_vn_free, color="r", scale=100)
quiver(gnss_lon, gnss_lat, gnss_pred_ve_fix, gnss_pred_vn_fix, color="b", scale=100)
axis("equal")
show()


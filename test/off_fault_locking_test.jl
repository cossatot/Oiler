using Pkg; Pkg.activate("..")

using Revise

using PyPlot

using CSV
using JSON
using DataFrames

using Oiler



block_df = Oiler.IO.gis_vec_file_to_df("./test_data/locking_test_blocks.geojson")


vel_df = Oiler.IO.gis_vec_file_to_df("./test_data/locking_test_points.csv")
vel_df[!, :fix] .= "1111"
vels = [Oiler.IO.vel_from_row(vel_df[i,:]; vel_type="GNSS") for i in 1:size(vel_df, 1)]#[1:10:end]

pole_sphere = Oiler.PoleSphere(lon=-25, lat=0.0, rotrate=-0.1,
                        fix="E", mov="W")

pole_cart = Oiler.pole_sphere_to_cart(pole_sphere)
ref_pole = Oiler.PoleCart(x=0.0, y=0.0, z=0.0, fix="1111", mov="E")

#poles = [pole_cart, ref_pole]

fault_df = Oiler.IO.load_faults_from_gis_files("./test_data/locking_test_faults.geojson")
fault_df[!,:v_rl] .= 0.0
fault_df[!,:e_rl] .= 0.0
fault_df[!,:v_ex] .= 0.0
fault_df[!,:e_ex] .= 0.0
fault_df[!,:usd] .= 0.0
fault_df[!,:lsd] .= 20.0
fault_df[!,:dip] .= 30.0

faults = Oiler.IO.process_faults_from_df(fault_df, )

full_fault = filter(x->x.fid=="full", faults)
ff_vels = Oiler.IO.make_vels_from_faults(full_fault)
segmented_faults = filter(x->x.fid!="full", faults)
sf_vels = Oiler.IO.make_vels_from_faults(segmented_faults)
alt_faults = segmented_faults[1:2:end]
olt_faults = segmented_faults[2:2:end]
af_vels = sf_vels[1:2:end]
non_fault_bounds = Oiler.IO.get_non_fault_block_bounds(block_df, alt_faults)
boundary_faults = map(Oiler.Boundaries.boundary_to_fault, non_fault_bounds)
bound_vels = Oiler.IO.make_vels_from_faults(boundary_faults)

vel_groups = Oiler.group_vels_by_fix_mov(vels)
#vel_groups[("E","W")] = Array{Oiler.VelocityVectorSphere}[]

ff_vel_groups = Oiler.group_vels_by_fix_mov(vcat(vels, ff_vels[1]))
sf_vel_groups = Oiler.group_vels_by_fix_mov(vcat(vels, sf_vels[1]))
af_vel_groups = Oiler.group_vels_by_fix_mov(vcat(vels, af_vels[1], bound_vels))


function reverse_fault(fault)
    new_dip_dir = Oiler.Faults.get_closest_dir(
        Oiler.Faults.get_trace_dip_trend_rhr(reverse(fault.trace, dims=1))
    )
    Oiler.Fault(
        trace=reverse(fault.trace, dims=1),
        dip=fault.dip,
        dip_dir = new_dip_dir,
        extension_rate=fault.extension_rate,
        extension_err=fault.extension_err,
        dextral_rate=fault.dextral_rate,
        cde=fault.cde,
        lsd=fault.lsd,
        usd=fault.usd,
        name=fault.name,
        hw=fault.fw,
        fw=fault.hw,
        fid=fault.fid,
        check_trace=false,
    )
end

bound_faults = map(reverse_fault, boundary_faults)



vd_segmented = Oiler.Solver.set_up_block_inv_no_constraints(
    sf_vel_groups,
    faults=segmented_faults,
)


vd_alternating = Oiler.Solver.set_up_block_inv_no_constraints(
    af_vel_groups,
    faults=alt_faults,
    bound_faults=map(reverse_fault, olt_faults),
)

vd_full = Oiler.Solver.set_up_block_inv_no_constraints(
    sf_vel_groups,
    faults=full_fault,
)


poles = Dict(
    ("1111", "E")=>ref_pole,
    ("1111", "W")=>pole_cart,
    ("E","W")=>pole_cart,
)

full_results = Dict(
    "predicted_vels"=>Oiler.ResultsAnalysis.predict_model_velocities(ff_vel_groups, vd_full, poles),
)
#Oiler.Plots.plot_results_map(full_results, sf_vel_groups, full_fault)
#ylim([1., 3.])
#title("full results")

seg_results = Dict(
    "predicted_vels"=>Oiler.ResultsAnalysis.predict_model_velocities(sf_vel_groups, vd_segmented, poles),
)

Oiler.Plots.plot_results_map(seg_results, sf_vel_groups, alt_faults)
for fault in boundary_faults
    plot(fault.trace[:, 1], fault.trace[:, 2], "b-", lw=0.5)
end
ylim([1., 3.])
title("all fault segments")


alt_results = Dict(
    "predicted_vels"=>Oiler.ResultsAnalysis.predict_model_velocities(af_vel_groups, vd_alternating, poles),
)

Oiler.Plots.plot_results_map(alt_results, af_vel_groups, alt_faults)
#for fault in boundary_faults
#    plot(fault.trace[:, 1], fault.trace[:, 2], "b-", lw=0.5)
#end
for fault in bound_faults
    plot(fault.trace[:, 1], fault.trace[:, 2], color="orange", ls=":", lw=1.5)
end
ylim([1., 3.])
title("alt faults and bound faults")

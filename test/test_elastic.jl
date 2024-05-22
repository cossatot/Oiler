using Test

using Oiler

rtol = 0.01



function test_okada_displacements_per_fault()

    ss = 1.0
    ds = 1.0
    ts = 1.0

    fault = Oiler.Fault(trace = [-81.905471 12.945880; -81.554797 12.114126],
        dip = 20.0, dip_dir = "S", lsd = 20.0)

    lons = [-81.739]
    lats = [12.527]

    n_gnss = length(lons)
    xp, yp = fault_oblique_merc(fault, lons, lats)
    xg, yg = xp[1:n_gnss], yp[1:n_gnss] # gnss
    sx1, sy1, sx2, sy2 = xp[end-1], yp[end-1], xp[end], yp[end] # fault 

    # format okada
    d = fault_to_okada(fault, sx1, sy1, sx2, sy2)

    # calc okada partials
    elastic_floor = 0.0
    offset_comps = okada(d, ss, ts, ds, xg, yg; floor = elastic_floor)

    x = offset_comps[1] + offset_comps[4] + offset_comps[7]
    y = offset_comps[2] + offset_comps[5] + offset_comps[8]
    z = offset_comps[3] + offset_comps[6] + offset_comps[9]

    x, y, z
end


function test_calc_locking_effects_per_fault_1()
    ff = Oiler.Fault(trace = [-81.905471 12.945880; -81.554797 12.114126],
        dip = 20.0, dip_dir = "S", lsd = 20.0)

    lons = [-81.739]
    lats = [12.527]

    pp = Oiler.Elastic.calc_okada_locking_effects_per_fault(ff, lons, lats;
        elastic_floor = 0.0)[1]

    # no good way of verifying this right now; doesn't match Meade and Loveless
    # because of the oblique mercator projections are different
    # and perhaps other reasons
    @test isapprox(pp, [1.11344e8 1.0689e9 4.68748e9
            -5.3601e9 -8.55279e8 -3.39197e8
            0.0 0.0 0.0]; rtol = rtol)
end


function test_calc_locking_effects_segmented_fault_1()

    # Should match test_calc_lockin_effects_per_fault_1()

    ff = Oiler.Fault(trace=[-81.905471 12.945880; 
                            -81.73 12.53;
                            -81.554797 12.114126],
                     dip=20., dip_dir="S", lsd=20.)

    lons = [-81.739]
    lats = [12.527]

    pp = Oiler.Elastic.calc_locking_effects_segmented_fault(ff, lons, lats; 
            elastic_floor=0.)[1]

    # no good way of verifying this right now; doesn't match Meade and Loveless
    # because of the oblique mercator projections are different
    # and perhaps other reasons
    @test isapprox(pp, [1.11344e8   1.0689e9    4.68748e9;
                        -5.3601e9   -8.55279e8  -3.39197e8;
                         0.0         0.0         0.0]; rtol=rtol)
end


function test_calc_locking_effects_1()
    ff = Oiler.Fault(trace=[-81.905471 12.945880; -81.554797 12.114126],
                     dip=20., dip_dir="S", lsd=20., hw="a", fw="b")
    vv = Oiler.VelocityVectorSphere(lon=-81.739, lat=12.527, vel_type="GNSS",
                                    ve=0., vn=0., fix="a", mov="b")

    fv = Oiler.Faults.fault_to_vel(ff)

    vels = convert(Array{VelocityVectorSphere}, [fv, vv])

    vel_groups = Oiler.IO.group_vels_by_fix_mov(vels)

    le = Oiler.Elastic.calc_locking_effects([ff], vel_groups; elastic_floor=0.)

    @test isapprox(le[4:6, 1:3], [1.11344e8   1.0689e9    4.68748e9;
                        -5.3601e9   -8.55279e8  -3.39197e8;
                         0.0         0.0         0.0]; rtol=rtol)

end


function test_get_fault_tris()
    new_thrust_coords = [-121.0 37.4
        -121.5 37.2
        -121.5 36.6
        -121.0 36.2]

    new_thrust = Oiler.Fault(trace = new_thrust_coords, dip = 15.0, lsd = 10.0, dip_dir = "E")

    tris = Oiler.Elastic.get_fault_tris(new_thrust)

    fault_tris = [Oiler.Tris.Tri([-121.0, 36.2, 0.0], [-121.5, 36.6, 0.0], [-120.60571152819729, 36.0938117593902, -10.0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "", "", "")
        Oiler.Tris.Tri([-121.10368477385849, 36.493802239447874, -10.0], [-120.60571152819729, 36.0938117593902, -10.0], [-121.5, 36.6, 0.0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "", "", "")
        Oiler.Tris.Tri([-121.5, 36.6, 0.0], [-121.5, 37.2, 0.0], [-121.10368477385849, 36.493802239447874, -10.0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "", "", "")
        Oiler.Tris.Tri([-121.100568518915, 37.09378777341805, -10.0], [-121.10368477385849, 36.493802239447874, -10.0], [-121.5, 37.2, 0.0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "", "", "")
        Oiler.Tris.Tri([-121.5, 37.2, 0.0], [-121.0, 37.4, 0.0], [-121.100568518915, 37.09378777341805, -10.0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "", "", "")
        Oiler.Tris.Tri([-120.599509046255, 37.29378290052262, -10.0], [-121.100568518915, 37.09378777341805, -10.0], [-121.0, 37.4, 0.0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "", "", "")]

    for (i, tri) in enumerate(tris)
        @test tri == fault_tris[i]
    end

end



@testset "test elastic.jl" begin

    test_calc_locking_effects_per_fault_1()
    test_calc_locking_effects_segmented_fault_1()
    test_calc_locking_effects_1()

end
using Test

using Oiler

rtol = 0.01



function test_calc_locking_effects_per_fault_1()
    ff = Oiler.Fault(trace=[-81.905471 12.945880; -81.554797 12.114126],
                     dip=20., dip_dir="S", lsd=20.)

    lons = [-81.739]
    lats = [12.527]

    pp = Oiler.Elastic.calc_locking_effects_per_fault(ff, lons, lats)[1]

    # no good way of verifying this right now; doesn't match Meade and Loveless
    # because of the oblique mercator projections are different
    # and perhaps other reasons
    @test isapprox(pp, [1.11344e8   1.0689e9    4.68748e9;
                        -5.3601e9   -8.55279e8  -3.39197e8;
                         0.0         0.0         0.0]; rtol=rtol)
end


function test_calc_locking_effects_segmented_fault_1()

    # Should match test_calc_lockin_effects_per_fault_1()

    ff = Oiler.Fault(trace=[-81.905471 12.945880; 
                            -81.73 12.53;
                            -81.554797 12.114126],
                     dip=20., dip_dir="S", lsd=20.)

    lons = [-81.739]
    lats = [12.527]

    pp = Oiler.Elastic.calc_locking_effects_segmented_fault(ff, lons, lats)[1]

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

    le = Oiler.Elastic.calc_locking_effects([ff], vel_groups)

    @test isapprox(le[4:6, 1:3], [1.11344e8   1.0689e9    4.68748e9;
                        -5.3601e9   -8.55279e8  -3.39197e8;
                         0.0         0.0         0.0]; rtol=rtol)

end


@testset "test elastic.jl" begin

    test_calc_locking_effects_per_fault_1()
    test_calc_locking_effects_segmented_fault_1()
    test_calc_locking_effects_1()

end
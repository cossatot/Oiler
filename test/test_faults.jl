using Test

using Oiler


function test_fault_slip_rate_to_ve_vn()
    sqrt2 = sqrt(2.0)
    test_dict = Dict((1.0, 0.0, 0.0) => (0.0, 1.0),
        (1.0, 0.0, 45.0) => (sqrt2 / 2.0, sqrt2 / 2.0),
        (1.0, 0.0, 90.0) => (1.0, 0.0),
        (1.0, 0.0, 135.0) => (sqrt2 / 2.0, -sqrt2 / 2.0), (0.0, 1.0, 0.0) => (-1.0, 0.0),
        (0.0, 1.0, 45.0) => (-sqrt2 / 2.0, sqrt2 / 2),
        (0.0, 1.0, 90.0) => (0.0, 1.0),
        (0.0, 1.0, 135.0) => (sqrt2 / 2.0, sqrt2 / 2.0), (1.0, 1.0, 0.0) => (-1.0, 1),
        (-1.0, 1.0, 0.0) => (-1, -1),
        (1.0, -1.0, 0.0) => (1, 1),
        (-1.0, -1.0, 0.0) => (1, -1), (1.0, 1.0, 45.0) => (0.0, sqrt2),
        (-1.0, 1.0, 45.0) => (-sqrt2, 0.0),
        (1.0, -1.0, 45.0) => (sqrt2, 0.0),
        (-1.0, -1.0, 45.0) => (0.0, -sqrt2),
        (1.0, 1.0, 90.0) => (1.0, 1.0),
        (-1.0, 1.0, 90.0) => (-1.0, 1.0),
        (1.0, -1.0, 90.0) => (1.0, -1.0),
        (-1.0, -1.0, 90.0) => (-1.0, -1.0),
        (1.0, 1.0, 135.0) => (sqrt2, 0.0),
        (-1.0, 1.0, 135.0) => (0.0, sqrt2),
        (1.0, -1.0, 135.0) => (0.0, -sqrt2),
        (-1.0, -1.0, 135.0) => (-sqrt2, 0.0))

    for (inputs, outputs) in test_dict
        @test isapprox(collect(Oiler.fault_slip_rate_to_ve_vn(inputs...)),
            collect(outputs))
    end
end


function test_ve_vn_to_fault_slip_rate()
    sqrt2 = sqrt(2.0)
    test_dict = Dict(
        (0.0, 1.0, 0.0) => (1.0, 0.0),
        (sqrt2 / 2.0, sqrt2 / 2.0, 45.0) => (1.0, 0.0),
        (1.0, 0.0, 90.0) => (1.0, 0.0),
        (sqrt2 / 2.0, -sqrt2 / 2.0, 135.0) => (1.0, 0.0), (-1.0, 0.0, 0.0) => (0.0, 1.0),
        (-sqrt2 / 2.0, sqrt2 / 2, 45.0) => (0.0, 1.0),
        (0.0, 1.0, 90.0) => (0.0, 1.0),
        (0.0, 1.0, 135.0) => (sqrt2 / 2.0, sqrt2 / 2.0), (1.0, 1.0, 0.0) => (-1.0, 1),
        (-1.0, 1.0, 0.0) => (-1, -1),
        (1.0, -1.0, 0.0) => (1, 1),
        (-1.0, -1.0, 0.0) => (1, -1), (0.0, sqrt2, 45.0) => (1.0, 1.0),
        (-sqrt2, 0.0, 45.0) => (-1.0, 1.0),
        (sqrt2, 0.0, 45.0) => (1.0, -1.0),
        (0.0, -sqrt2, 45.0) => (-1.0, -1.0),
        (1.0, 1.0, 90.0) => (1.0, 1.0),
        (-1.0, 1.0, 90.0) => (-1.0, 1.0),
        (1.0, -1.0, 90.0) => (1.0, -1.0),
        (-1.0, -1.0, 90.0) => (-1.0, -1.0),
        (sqrt2, 0.0, 135.0) => (1.0, 1.0),
        (0.0, sqrt2, 135.0) => (-1.0, 1.0),
        (0.0, -sqrt2, 135.0) => (1.0, -1.0),
        (-sqrt2, 0.0, 135.0) => (-1.0, -1.0)
    )

    for (inputs, outputs) in test_dict
        @test isapprox(collect(Oiler.Faults.ve_vn_to_fault_slip_rate(inputs...)),
            collect(outputs))
    end
end

quien_sabe_coords = [-121.37575 36.93768
    -121.33534 36.89357
    -121.29735 36.85596
    -121.25795 36.81047
    -121.2073 36.75499];


function test_get_midpoint_old()
    @test Oiler.Faults.get_midpoint_old(quien_sabe_coords) == (-121.29735, 36.85596)
end


function test_get_midpoint()
    @test Oiler.Faults.get_midpoint(quien_sabe_coords) == [-121.28985450382699 36.84731110410578]
end


function test_fault_to_vel_1()
    qs_1 = Oiler.Fault(trace = quien_sabe_coords, dip = 89.0, dip_dir = "NE",
        name = "Quien Sabe NE dip", hw = "gab", fw = "hol", dextral_rate = 35.0)

    vel = Oiler.Faults.fault_to_vel(qs_1)

    @test vel.lon == -121.28985450382699
    @test vel.lat == 36.84731110410578
    @test vel.ve == -20.77427972024876
    @test vel.vn == 28.167877131670057
    @test vel.vu == 0.0
    @test vel.ee == 0.0
    @test vel.en == 0.0
    @test vel.eu == 0.0
    @test vel.fix == "gab"
    @test vel.mov == "hol"
    @test vel.name == "Quien Sabe NE dip"
    @test vel.vel_type == "fault"
end


function test_fault_to_vel_2()
    qs_2 = Oiler.Fault(trace = quien_sabe_coords, dip = 89.0, dip_dir = "SW",
        name = "Quien Sabe S dip", hw = "hol", fw = "gab", dextral_rate = 35.0)

    vel = Oiler.Faults.fault_to_vel(qs_2)

    @test vel.lon == -121.28985450382699
    @test vel.lat == 36.84731110410578
    @test vel.ve == 20.786879116086872
    @test vel.vn == -28.158580514883763
    @test vel.vu == 0.0
    @test vel.ee == 0.0
    @test vel.en == 0.0
    @test vel.eu == 0.0
    @test vel.fix == "hol"
    @test vel.mov == "gab"
    @test vel.name == "Quien Sabe S dip"
    @test vel.vel_type == "fault"
end


function test_fault_to_vel_3()
    qs_3 = Oiler.Fault(trace = quien_sabe_coords, dip = 90.0, dip_dir = "V",
        name = "Quien Sabe V dip", hw = "hol", fw = "gab", dextral_rate = 35.0)

    vel = Oiler.Faults.fault_to_vel(qs_3)

    @test vel.lon == -121.28985450382699
    @test vel.lat == 36.84731110410578
    @test vel.ve == 20.786879116086872
    @test vel.vn == -28.158580514883763
    @test vel.vu == 0.0
    @test vel.ee == 0.0
    @test vel.en == 0.0
    @test vel.eu == 0.0
    @test vel.fix == "hol"
    @test vel.mov == "gab"
    @test vel.name == "Quien Sabe V dip"
    @test vel.vel_type == "fault"
end


function test_check_right_hand_rule_no_change()
    trace_in = quien_sabe_coords
    dip_dir = "SW"
    trace_out = Oiler.Faults.check_right_hand_rule(trace_in, dip_dir)
    @test trace_in == trace_out
end


function test_check_right_hand_rule_flip()
    trace_in = quien_sabe_coords
    dip_dir = "NE"
    trace_out = Oiler.Faults.check_right_hand_rule(trace_in, dip_dir)
    @test reverse(trace_in, dims = 1) == trace_out
end


function test_get_fault_slip_rate_from_pole_reverse_no_change()
    fault = Oiler.Fault(trace = [5.0 0.0; -5.0 0.0], dip_dir = "N", dip = 30.0, hw = "N", fw = "S")
    pole = Oiler.PoleCart(x = -4.312474367375016e-10,
        y = -1.5696123057604772e-10, z = 0.0, mov = "S", fix = "N")

    fsr = Oiler.Faults.get_fault_slip_rate_from_pole(fault, pole)

    @test isapprox(fsr[1], 0.0, atol = 1e-10)
    @test isapprox(fsr[2], -1.0, atol = 1e-10)
end


function test_get_fault_slip_rate_from_pole_reverse_flip_pole()
    fault = Oiler.Fault(trace = [-5.0 0.0; 5.0 0.0], dip_dir = "S", dip = 30.0, hw = "S", fw = "N")
    pole = Oiler.PoleCart(x = -4.312474367375016e-10,
        y = -1.5696123057604772e-10, z = 0.0, mov = "S", fix = "N")

    fsr = Oiler.Faults.get_fault_slip_rate_from_pole(fault, pole)

    @test isapprox(fsr[1], 0.0, atol = 1e-10)
    @test isapprox(fsr[2], -1.0, atol = 1e-10)
end


function test_get_fault_slip_rate_from_pole_mismatch()
    fault = Oiler.Fault(trace = [-5.0 0.0; 5.0 0.0], dip_dir = "S", dip = 30.0, hw = "N", fw = "S")
    pole = Oiler.PoleCart(x = -4.312474367375016e-10,
        y = -1.5696123057604772e-10, z = 0.0, mov = "S", fix = "A")

    fsr = Oiler.Faults.get_fault_slip_rate_from_pole(fault, pole)

    @test isapprox(fsr[1], 0.0, atol = 1e-10)
    @test isapprox(fsr[2], 1.0, atol = 1e-10)
end



function test_get_fault_slip_rate_from_pole_strikeslip_no_change()
    fault = Oiler.Fault(trace = [0.0 -5.0; 0.0 5.0], dip_dir = "E", dip = 89.0, hw = "E", fw = "W")
    pole = Oiler.PoleCart(x = -4.312474367375016e-10,
        y = -1.5696123057604772e-10, z = 0.0, mov = "W", fix = "E")

    fsr = Oiler.Faults.get_fault_slip_rate_from_pole(fault, pole)

    @test isapprox(fsr[1], 1.0, atol = 1e-10)
    @test isapprox(fsr[2], 0.0, atol = 1e-10)
end


function test_get_fault_slip_rate_from_pole_strikeslip_flip()
    fault = Oiler.Fault(trace = [0.0 -5.0; 0.0 5.0], dip_dir = "W", dip = 89.0, hw = "W", fw = "E")
    pole = Oiler.PoleCart(x = -4.312474367375016e-10,
        y = -1.5696123057604772e-10, z = 0.0, mov = "W", fix = "E")

    fsr = Oiler.Faults.get_fault_slip_rate_from_pole(fault, pole)

    @test isapprox(fsr[1], 1.0, atol = 1e-10)
    @test isapprox(fsr[2], 0.0, atol = 1e-10)
end


function test_fault_to_vels() end


function test_fault_oblique_merc()

    ff = Oiler.Fault(trace = [83.0 32.0; 83.5 32.8], dip = 60.0, dip_dir = "SE")

    lons = [100.0; 101.0; 122.0]
    lats = [20.0; 20.0; 80.5]

    x, y = Oiler.Faults.fault_oblique_merc(ff, lons, lats)

    # these values from faultobliquemerc.m by Meade and Loveless
    x_ = [-0.995670130919819; -0.985654810206086; -0.099467883321707]
    y_ = [-0.342577075025567; -0.356822693071572; 0.271564430987158]

    @test isapprox(x[1:3], x_ .* Oiler.EARTH_RAD_KM)
    @test isapprox(y[1:3], y_ .* Oiler.EARTH_RAD_KM)
end


function test_strike_rot_matrix_xd()
    xd_partials = [0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0]
    strike = 45.0

    srt = Oiler.Faults.build_strike_rot_matrix(strike)

    v = srt * xd_partials * [0.0; 1.0; 0.0]
    @test isapprox(v, [-0.7071067811865475; 0.7071067811865475; 0.0])
end


function test_build_Pf_vert() end


function test_build_Pf_dip() end




function test_partials_and_rotation_matrices_1()
    # not sure what is supposed to happen in this function

    fault = Oiler.Fault(trace = [-5.0 0.0; 5.0 0.0], dip_dir = "S", dip = 30.0)

    pp = Oiler.Elastic.calc_locking_effects_per_fault(fault, [0.0], [-0.0001])[1]

    pole = [-4.312474367375016e-10; -1.5696123057604772e-10; 0.0]

    pp * pole
end


function test_build_velocity_projection_matrix_1()
    Pf = Oiler.build_velocity_projection_matrix(90.0, 90.0)

    @test isapprox(Pf * [1.0; 0.0; 0.0], [1.0; 0.0; 0.0])
end


function test_project_fault_trace()
    new_thrust_coords = [-121.0 37.4
        -121.5 37.2
        -121.5 36.6
        -121.0 36.2]

    new_thrust = Oiler.Fault(trace = new_thrust_coords, dip = 15.0, lsd = 10.0, dip_dir = "E")

    new_trace = Oiler.Faults.project_fault_trace(new_thrust)

    new_trace_ans = [-120.5840782861774 36.20045428065595; 
    -121.08193195903539 36.6004436871955; 
    -121.07863183830526 37.200427589088946; 
    -120.57750984903879 37.40042216620642]

    @test new_trace == new_trace_ans

end


@testset "test faults.jl" begin
    test_fault_to_vel_1()
    test_fault_to_vel_2()
    test_fault_to_vel_3()
    test_get_midpoint()
    test_get_midpoint_old()
    test_check_right_hand_rule_no_change()
    test_check_right_hand_rule_flip()
    test_get_fault_slip_rate_from_pole_reverse_no_change()
    test_get_fault_slip_rate_from_pole_reverse_flip_pole()
    test_get_fault_slip_rate_from_pole_mismatch()
    test_get_fault_slip_rate_from_pole_strikeslip_no_change()
    test_get_fault_slip_rate_from_pole_strikeslip_flip()
    test_fault_slip_rate_to_ve_vn()
    test_strike_rot_matrix_xd()
    test_build_velocity_projection_matrix_1()
    test_project_fault_trace()
end
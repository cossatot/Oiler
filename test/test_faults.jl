using Revise

using Test


using Oiler

quien_sabe_coords = [-121.37575 36.93768;
                     -121.33534 36.89357;
                     -121.29735 36.85596;
                     -121.25795 36.81047;
                     -121.2073 36.75499]

function test_fault_slip_rate_to_ve_vn()
    sqrt2 = sqrt(2.)
    test_dict = Dict((1., 0., 0.) => (0., 1.),
        (1., 0., 45.) => (sqrt2 / 2., sqrt2 / 2.),
        (1., 0., 90.) => (1., 0.),
        (1., 0., 135.) => (sqrt2 / 2., -sqrt2 / 2.),
        
        (0., 1., 0.) => (-1., 0.),
        (0., 1., 45.) => (-sqrt2 / 2., sqrt2 / 2),
        (0., 1., 90.) => (0., 1.),
        (0., 1., 135.) => (sqrt2 / 2., sqrt2 / 2.),

        (1., 1., 0.) => (-1., 1),
        (-1., 1., 0.) => (-1, -1),
        (1., -1., 0.) => (1, 1),
        (-1., -1., 0.) => (1, -1),

        (1., 1., 45.) => (0., sqrt2),
        (-1., 1., 45.) => (-sqrt2, 0.),
        (1., -1., 45.) => (sqrt2, 0.),
        (-1., -1., 45.) => (0., -sqrt2),
        (1., 1., 90.) => (1., 1.),
        (-1., 1., 90.) => (-1., 1.),
        (1., -1., 90.) => (1., -1.),
        (-1., -1., 90.) => (-1., -1.),
        (1., 1., 135.) => (sqrt2, 0.),
        (-1., 1., 135.) => (0., sqrt2),
        (1., -1., 135.) => (0., -sqrt2),
        (-1., -1., 135.) => (-sqrt2, 0.))

    for (inputs, outputs) in test_dict
        @test isapprox(collect(Oiler.fault_slip_rate_to_ve_vn(inputs...)), 
                       collect(outputs))

    end
end

test_fault_slip_rate_to_ve_vn()

# @test Oiler.get_midpoint(quien_sabe_coords) == (-121.29735, 36.85596)

qs_1 = Oiler.Fault(trace = quien_sabe_coords, dip = 89., dip_dir = "NE", 
    name = "Quien Sabe NE dip", hw = "gab", fw = "hol", dextral_rate = 35.);
qs_2 = Oiler.Fault(trace = quien_sabe_coords, dip = 89., dip_dir = "S",
    name = "Quien Sabe S dip", hw = "hol", fw = "gab", dextral_rate = 35.);
qs_3 = Oiler.Fault(trace = quien_sabe_coords, dip = 90., dip_dir = "V", 
    name = "Quien Sabe V dip", hw = "hol", fw = "gab", dextral_rate = 35.);

Oiler.fault_to_vel(qs_1);


function test_fault_oblique_merc()

    ff = Oiler.Fault(trace = [83. 32.; 83.5 32.8], dip = 60., dip_dir = "SE")

    lons = [100.; 101.; 122.]
    lats = [20.; 20.; 80.5]

    x, y = Oiler.Faults.fault_oblique_merc(ff, lons, lats)

    # these values from faultobliquemerc.m by Meade and Loveless
    x_ = [-0.995670130919819; -0.985654810206086; -0.099467883321707]
    y_ = [-0.342577075025567; -0.356822693071572;  0.271564430987158]

    @test isapprox(x[1:3], x_ .* Oiler.EARTH_RAD_KM)
    @test isapprox(y[1:3], y_ .* Oiler.EARTH_RAD_KM)
end

#test_fault_oblique_merc()


function test_strike_rot_matrix_xd()
    xd_partials = [0. 0. 0.; 0. 1. 0.; 0. 0. 0.]
    strike = 45.

    srt = Oiler.Faults.build_strike_rot_matrix(strike)

    v = srt * xd_partials * [0.; 1.; 0.]
    @test isapprox(v, [-0.7071067811865475; 0.7071067811865475; 0.])
end

test_strike_rot_matrix_xd()

function test_partials_and_rotation_matrices_1()
    fault = Oiler.Fault(trace = [-5. 0.; 5. 0.], dip_dir = "S", dip = 30.)
    
    pp = Oiler.Elastic.calc_locking_effects_per_fault(fault, [0.], [-0.0001])[1]

    pole = [-4.312474367375016e-10; -1.5696123057604772e-10; 0.]

    pp * pole
end

#test_partials_and_rotation_matrices_1()

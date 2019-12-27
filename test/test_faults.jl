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
    name = "Quien Sabe NE dip", hw = "gab", fw = "hol", dextral_rate = 35.)
qs_2 = Oiler.Fault(trace = quien_sabe_coords, dip = 89., dip_dir = "S",
    name = "Quien Sabe S dip", hw = "hol", fw = "gab", dextral_rate = 35.)
qs_3 = Oiler.Fault(trace = quien_sabe_coords, dip = 90., dip_dir = "V", 
    name = "Quien Sabe V dip", hw = "hol", fw = "gab", dextral_rate = 35.)

Oiler.fault_to_vel(qs_1)
using Test

using Oiler
using Oiler: VelocityVectorSphere

using Oiler.Utils: get_gnss_vels

vel_groups = Dict(("a", "b") => [
    VelocityVectorSphere(lon = 0.0,
                         lat = 0.05,
                         ve = 0.5,
                         vn = 0.,
                         fix = "a",
                         mov = "b",
                         name = "f1",
                         vel_type = "fault"),
    VelocityVectorSphere(lon = 0.05,
                         lat = 0.,
                         ve = 0.53,
                         vn = 0.,
                         fix = "a",
                         mov = "b",
                         name = "f2",
                         vel_type = "fault")
    ],

    ("b", "c") => [
    VelocityVectorSphere(lon = 0.7,
                         lat = 0.05,
                         ve = 1.5,
                         vn = 0.,
                         fix = "c",
                         mov = "b",
                         name = "f3",
                         vel_type = "fault")
],
    ("r", "a") => [
    VelocityVectorSphere(lon = -1.,
                         lat = 0.05,
                         ve = 0.0,
                         vn = 0.,
                         fix = "r",
                         mov = "a",
                         name = "g3",
                         vel_type = "GNSS"),

    VelocityVectorSphere(lon = 0.05,
                         lat = 0.,
                         ve = 0.15,
                         vn = 0.,
                         fix = "r",
                         mov = "a",
                         name = "g4",
                         vel_type = "GNSS")
    ],
    ("r", "c") => [
    VelocityVectorSphere(lon = 0.77,
                         lat = 0.05,
                         ve = 1.9,
                         vn = 0.,
                         fix = "r",
                         mov = "c",
                         name = "g5",
                         vel_type = "GNSS")
])


faults = [Oiler.Fault(trace=[0. 0.; 1. 1.], dip=89., dip_dir="SW", hw="a",
                    fw="b", name="ff1"),
          Oiler.Fault(trace=[5. 5.; 1. 1.], dip=89., dip_dir="SW", hw="a",
                    fw="b", name="ff2"),
          Oiler.Fault(trace=[5. 5.; 1. 1.], dip=89., dip_dir="NE", hw="c",
                    fw="b", name="ff3"),
          Oiler.Fault(trace=[5. 5.; 1. 1.], dip=89., dip_dir="SW", hw="b",
                    fw="c", name="ff4"),
]


function test_get_gnss_vels()
    gnss_vel_answer = [Dict("idx" => UnitRange{Int64}[10:12, 7:9],"vel" => VelocityVectorSphere(lon = -1.,
    lat = 0.05,
    ve = 0.0,
    vn = 0.,
    fix = "r",
    mov = "a",
    name = "g3",
    vel_type = "GNSS"))
    Dict("idx" => UnitRange{Int64}[13:15, 7:9],"vel" => VelocityVectorSphere(lon = 0.05,
    lat = 0.,
    ve = 0.15,
    vn = 0.,
    fix = "r",
    mov = "a",
    name = "g4",
    vel_type = "GNSS"))
    Dict("idx" => UnitRange{Int64}[16:18, 10:12],"vel" => VelocityVectorSphere(lon = 0.77,
    lat = 0.05,
    ve = 1.9,
    vn = 0.,
    fix = "r",
    mov = "c",
    name = "g5",
    vel_type = "GNSS"))]
    
    gnss_vels = get_gnss_vels(vel_groups)

    @test gnss_vels == gnss_vel_answer
end

test_get_gnss_vels()

function test_get_fault_vels()
    get_fault_vels_answer = [
    Dict{Any,Any}("idx" => UnitRange{Int64}[1:3, 1:3],"fault" => VelocityVectorSphere(
    lon=0.0,
    lat=0.05,
    ve= 0.5,
    vn= 0.0,
    vu= 0.0,
    ee= 0.0,
    en= 0.0,
    eu= 0.0,
    fix= "a",
    mov= "b",
    name= "f1",
    vel_type= "fault")
    ),
    Dict{Any,Any}("idx" => UnitRange{Int64}[4:6, 1:3],"fault" => VelocityVectorSphere(
    lon=0.05,
    lat=0.0,
    ve= 0.53,
    vn= 0.0,
    vu= 0.0,
    ee= 0.0,
    en= 0.0,
    eu= 0.0,
    fix= "a",
    mov= "b",
    name= "f2",
    vel_type= "fault")
    ),
    Dict{Any,Any}("idx" => UnitRange{Int64}[7:9, 4:6],"fault" => VelocityVectorSphere(
    lon=0.7,
    lat=0.05,
    ve= 1.5,
    vn= 0.0,
    vu= 0.0,
    ee= 0.0,
    en= 0.0,
    eu= 0.0,
    fix= "c",
    mov= "b",
    name= "f3",
    vel_type= "fault")
    ),
    ]

    fault_vels = Oiler.Utils.get_fault_vels(vel_groups)
    @test fault_vels == get_fault_vels_answer
end

test_get_fault_vels()


function test_group_faults()
    vg_keys = sort(collect(Tuple(keys(vel_groups))))
    fault_groups = Oiler.Utils.group_faults(faults, vg_keys)

    group_faults_answer = Dict(
        ("a","b") => [faults[1], faults[2]],
        ("b","c") => [faults[3], faults[4]]
    )

    @test fault_groups == group_faults_answer
end

test_group_faults()
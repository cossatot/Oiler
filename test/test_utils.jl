using Test

using Oiler
using Oiler: VelocityVectorSphere

using Oiler.Utils: get_gnss_vels

vel_groups = Dict(("a", "b") => [
    VelocityVectorSphere(lond = 0.0,
                         latd = 0.05,
                         ve = 0.5,
                         vn = 0.,
                         fix = "a",
                         mov = "b",
                         name = "f1",
                         vel_type = "fault"),
    VelocityVectorSphere(lond = 0.05,
                         latd = 0.,
                         ve = 0.53,
                         vn = 0.,
                         fix = "a",
                         mov = "b",
                         name = "f2",
                         vel_type = "fault")
    ],

    ("b", "c") => [
    VelocityVectorSphere(lond = 0.7,
                         latd = 0.05,
                         ve = 1.5,
                         vn = 0.,
                         fix = "a",
                         mov = "b",
                         name = "f3",
                         vel_type = "fault")
],
    ("r", "a") => [
    VelocityVectorSphere(lond = -1.,
                         latd = 0.05,
                         ve = 0.0,
                         vn = 0.,
                         fix = "a",
                         mov = "b",
                         name = "g3",
                         vel_type = "GNSS"),

    VelocityVectorSphere(lond = 0.05,
                         latd = 0.,
                         ve = 0.15,
                         vn = 0.,
                         fix = "a",
                         mov = "b",
                         name = "g4",
                         vel_type = "GNSS")
    ],
    ("r", "c") => [
    VelocityVectorSphere(lond = 0.77,
                         latd = 0.05,
                         ve = 1.9,
                         vn = 0.,
                         fix = "r",
                         mov = "c",
                         name = "g5",
                         vel_type = "GNSS")
])

function test_get_gnss_vels()
    gnss_vel_answer = [Dict("idx" => UnitRange{Int64}[1:3, 1:3],"vel" => VelocityVectorSphere(lond = -1.,
    latd = 0.05,
    ve = 0.0,
    vn = 0.,
    fix = "a",
    mov = "b",
    name = "g3",
    vel_type = "GNSS"))
    Dict("idx" => UnitRange{Int64}[4:6, 1:3],"vel" => VelocityVectorSphere(lond = 0.05,
    latd = 0.,
    ve = 0.15,
    vn = 0.,
    fix = "a",
    mov = "b",
    name = "g4",
    vel_type = "GNSS"))
    Dict("idx" => UnitRange{Int64}[10:12, 7:9],"vel" => VelocityVectorSphere(lond = 0.77,
    latd = 0.05,
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
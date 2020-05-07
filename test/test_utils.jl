using Test

using Oiler
using Oiler: VelocityVectorSphere

using Oiler.Utils: get_gnss_vels


function test_diagonalize_matrices()
end

function test_make_digraph_from_vels_many_args()
end

function test_make_digraph_from_vels_array()
end

function test_make_digraph_from_poles()
end

function test_make_digraph_from_tuples()
    keys = [("a", "b"), ("b", "c"), ("c", "a"), ("r", "a"), ("r", "c")]
    dg = Dict("c"=>["a"], "b"=>["c"], "r"=>["a", "c"], "a"=>["b"])
    @test Oiler.Utils.make_digraph_from_tuples(keys) == dg
end

test_make_digraph_from_tuples()

function test_make_ugraph_from_digraph()
    dg = Dict("c"=>["a"], "b"=>["c"], "r"=>["a", "c"], "a"=>["b"])
    ug = Dict(
        "c"=>["a", "b", "r"],
        "b"=>["c", "a"],
        "r"=>["a", "c"],
        "a"=>["c", "r", "b"]
        )
    @test Oiler.Utils.make_ugraph_from_digraph(dg) == ug
end

test_make_ugraph_from_digraph()

function test_find_tricycles()
    ug = Dict(
        "c"=>["a", "b", "r"],
        "b"=>["c", "a"],
        "r"=>["a", "c"],
        "a"=>["c", "r", "b"]
        )
    tris = [["c", "a", "r"], ["c", "a", "b"]]

    @test Oiler.Utils.find_tricycles(ug) == tris
end

test_find_tricycles()


function test_get_cycle_inds_fw()
    keys = [("a", "b"), ("b", "c"), ("c", "a"), ("r", "a"), ("r", "c")]
    cycle_tup = ("c", "a")
    cycle_ind = Dict("ind"=>3, "val"=>1.0)
    @test Oiler.Utils.get_cycle_inds(keys, cycle_tup) == cycle_ind
end

test_get_cycle_inds_fw()


function test_get_cycle_inds_rev()
    keys = [("a", "b"), ("b", "c"), ("c", "a"), ("r", "a"), ("r", "c")]
    cycle_tup = ("a", "c")
    cycle_ind = Dict("ind"=>3, "val"=>-1.0)
    @test Oiler.Utils.get_cycle_inds(keys, cycle_tup) == cycle_ind
end

test_get_cycle_inds_rev()


function test_find_vel_cycles()
    keys = [("a", "b"), ("b", "c"), ("c", "a"), ("r", "a"), ("r", "c")]
    vel_cycles = Dict(
        1 => Dict(
            ("c", "a") => Dict{String,Real}("ind"=>3,"val"=>1.0),
            ("r", "c") => Dict{String,Real}("ind"=>5,"val"=>1.0),
            ("a", "r") => Dict{String,Real}("ind"=>4,"val"=>-1.0)
        ),
        2 => Dict(
            ("b", "c") => Dict{String,Real}("ind"=>2,"val"=>1.0),
            ("c", "a") => Dict{String,Real}("ind"=>3,"val"=>1.0),
            ("a", "b") => Dict{String,Real}("ind"=>1,"val"=>1.0)
        )
    )
    
    @test Oiler.Utils.find_vel_cycles(keys) == vel_cycles
end

test_find_vel_cycles()


function test_flat()
end


function test_find_shortest_path()
end


function test_get_pole_path()
end

function test_get_path_euler_pole()
end


function test_predict_vels_from_poles_cart()
end


function test_predict_vels_from_poles_sphere()
end

function test_get_coords_from_vel_array()
end


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
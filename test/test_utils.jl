using Test

using SparseArrays

using Oiler
using Oiler: VelocityVectorSphere

using Oiler.Utils: get_gnss_vels, sparse_to_dict, dict_to_sparse


function test_diagonalize_matrices()
end


function test_sparse_to_dict()
    I = [1, 4, 3, 5]; J = [4, 7, 18, 9]; V = [1., 2., -5., 3.];
    S = sparse(I, J, V)
    sd = sparse_to_dict(S)

    sd_test = Dict((4, 7)  => 2.0,
                   (5, 9)  => 3.0,
                   (1, 4)  => 1.0,
                   (3, 18) => -5.0)

    @test sd == sd_test
end


function test_dict_to_sparse()
    I = [1, 4, 3, 5]; J = [4, 7, 18, 9]; V = [1., 2., -5., 3.];
    S_test = sparse(I, J, V)
    sd_test = Dict((4, 7)  => 2.0,
                   (5, 9)  => 3.0,
                   (1, 4)  => 1.0,
                   (3, 18) => -5.0)
    S = dict_to_sparse(sd_test)
    @test S == S_test
end


function test_make_digraph_from_vels_many_args()
end

function test_make_digraph_from_vels_array()
end

function test_make_digraph_from_vels_array_dup_fix_mov()
    na_1 = Oiler.VelocityVectorSphere( lon= -121.8, lat= 46.8, ve=16.6,
        vn= 0.25, vu= 0.0, ee= 1.0, en= 1.0, eu= 0.0,
        fix= "pa", mov= "na", name="na_1")

    jf_1 = Oiler.VelocityVectorSphere( lon= -127.4, lat= 45.3, ve=47.2,
        vn= -17.7, vu= 0.0, ee= 1.0, en= 1.0, eu= 0.0,
        fix= "pa", mov= "jf", name="jf_1")

    pa_pa = Oiler.VelocityVectorSphere( lon= -126.0, lat= 41.3, ve=43.4,
        vn= -19.3, vu= 0.0, ee= 1.0, en= 1.0, eu= 0.0,
        fix= "pa", mov= "pa", name="pa_pa")

    dg = Oiler.Utils.make_digraph_from_vels([na_1, jf_1, pa_pa])
    
  # @test_logs Oiler.Utils.make_digraph_from_vels([na_1, jf_1, pa_pa])
  @test dg["pa"] == ["na", "jf"]
end

function test_make_digraph_from_poles()
end

function test_make_digraph_from_tuples()
    keys = [("a", "b"), ("b", "c"), ("c", "a"), ("r", "a"), ("r", "c")]
    dg = Dict("c"=>["a"], "b"=>["c"], "r"=>["a", "c"], "a"=>["b"])
    @test Oiler.Utils.make_digraph_from_tuples(keys) == dg
end


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



function test_get_cycle_inds_fw()
    keys = [("a", "b"), ("b", "c"), ("c", "a"), ("r", "a"), ("r", "c")]
    cycle_tup = ("c", "a")
    cycle_ind = Dict("ind"=>3, "val"=>1.0)
    @test Oiler.Utils.get_cycle_inds(keys, cycle_tup) == cycle_ind
end



function test_get_cycle_inds_rev()
    keys = [("a", "b"), ("b", "c"), ("c", "a"), ("r", "a"), ("r", "c")]
    cycle_tup = ("a", "c")
    cycle_ind = Dict("ind"=>3, "val"=>-1.0)
    @test Oiler.Utils.get_cycle_inds(keys, cycle_tup) == cycle_ind
end



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



function test_group_faults()
    vg_keys = sort(collect(Tuple(keys(vel_groups))))
    fault_groups = Oiler.Utils.group_faults(faults, vg_keys)

    group_faults_answer = Dict(
        ("a","b") => [faults[1], faults[2]],
        ("b","c") => [faults[3], faults[4]]
    )

    @test fault_groups == group_faults_answer
end

function test_find_first_nz_per_row()

    A = sparse([0.0 2.3 0 0; 0. 0. 9. 22.; 2. 0 5. 7.])
    first_idxs = Oiler.Utils.find_first_nz_per_row(A)
    first_idxs_answer = [(1,2), (2,3), (3,1)]
    @test first_idxs == first_idxs_answer
end


function test_sort_sparse_matrix()
    A = sparse([0.0 2.3 0 0; 0. 0. 9. 22.; 2. 0 5. 7.])

    A_sort = Oiler.Utils.sort_sparse_matrix(A)
    A_sort_ans = sparse([2. 0. 5. 7.; 0. 2.3 0. 0.; 0. 0. 9. 22.])

    @test A_sort == A_sort_ans
end

function test_sort_segs()
    test_set = Set([
        [-1. -1.; 1. 1.],
        [1. 1.; 2. 2.],
        [2. 2.; 3. 3.],
        ]
    )

    sorted_ans = [-1. -1.; 1. 1.; 2. 2.; 3. 3.]

end

@testset "test utils.jl" begin

test_sparse_to_dict()
test_dict_to_sparse()
test_make_digraph_from_vels_array_dup_fix_mov()
test_make_digraph_from_tuples()
test_make_ugraph_from_digraph()
test_find_tricycles()
test_get_cycle_inds_fw()
test_get_cycle_inds_rev()
test_find_vel_cycles()
test_get_gnss_vels()
test_get_fault_vels()
test_group_faults()
end
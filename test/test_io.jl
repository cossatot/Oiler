using Test

using Oiler


# from https://stackoverflow.com/questions/66801702/
# deriving-equality-for-julia-structs-with-mutable-members
function struct_equals(a::S, b::S) where S
    for name in fieldnames(S)
        if getfield(a, name) != getfield(b, name)
            return false
        end
    end
    return true
end



test_gpkg = "./test_data/io_test.gpkg"
test_gj_blocks = "./test_data/io_test_blocks.geojson"

function test_gis_vec_file_to_df_gpkg_faults()
    fault_df = Oiler.IO.gis_vec_file_to_df(test_gpkg;
        layername="faults")

    @test size(fault_df) == (3, 15)
    @test fault_df[:, :fid] == [1, 2, 4]
    @test fault_df[:, :dip_dir] == ["E", "NW", "E"]
end

function test_gis_vec_file_to_df_gpkg_blocks()
    block_df = Oiler.IO.gis_vec_file_to_df(test_gpkg;
        layername="blocks")
end

function test_gis_vec_file_to_df_gpkg_gnss()
    gnss_df = Oiler.IO.gis_vec_file_to_df(test_gpkg;
        layername="gnss_vels")
end

function test_gis_vec_file_to_df_geojson_blocks()
    block_df = Oiler.IO.gis_vec_file_to_df(test_gpkg;
        layername="blocks")
    @test size(block_df, 1) == 4
end


function load_geodataframes()
    fault_df = Oiler.IO.gis_vec_file_to_df(test_gpkg; layername="faults")
    block_df = Oiler.IO.gis_vec_file_to_df(test_gpkg; layername="blocks")
    gnss_df = Oiler.IO.gis_vec_file_to_df(test_gpkg;layername="gnss_vels")

    (fault_df, block_df, gnss_df)
end


function test_row_to_fault()
    fault_df = Oiler.IO.gis_vec_file_to_df(test_gpkg; layername="faults")

    fault = Oiler.IO.row_to_fault(fault_df[1,:])

    fault_ans = Oiler.Fault(;trace=[0.49262599938633844 0.7071696809078425; 
                               0.6397098710932747 1.7367567828563972],
                               dip=60.,
                               dip_dir="E",
                               extension_rate=0.,
                               extension_err=1.,
                               dextral_rate=0.,
                               dextral_err=1.,
                               lsd=20.,
                               usd=0.,
                               name="f1",
                               hw="3",
                               fw="5",
                               fid=1
                               )

    @test struct_equals(fault, fault_ans)
end


function test_process_faults_from_df()

    fault_df = Oiler.IO.gis_vec_file_to_df(test_gpkg; layername="faults")

    faults = Oiler.IO.process_faults_from_df(fault_df)
    
    fault_ans_1 = Oiler.Fault(;trace=[0.49262599938633844 0.7071696809078425; 
                               0.6397098710932747 1.7367567828563972],
                               dip=60.,
                               dip_dir="E",
                               extension_rate=0.,
                               extension_err=1.,
                               dextral_rate=0.,
                               dextral_err=1.,
                               lsd=20.,
                               usd=0.,
                               name="f1",
                               hw="3",
                               fw="5",
                               fid=1
                               )

    fault_ans_2 = Oiler.Fault(;trace=[1.6490095424615623 2.233798832072941; 
                                    1.8163808447487657 1.4882357582481254; 
                                    0.9288057568620807 1.2904333100905214],
                               dip=30.,
                               dip_dir="NW",
                               extension_rate=1.,
                               extension_err=5.,
                               dextral_rate=3.,
                               dextral_err=2.,
                               lsd=15.,
                               usd=0.,
                               name="f2",
                               hw="1",
                               fw="2",
                               fid=2
                               )

    fault_ans_3 = Oiler.Fault(;trace=[0.9288057568620807 1.2904333100905214; 
                                    0.6397098710932747 1.7367567828563972],
                               dip=89.,
                               dip_dir="E",
                               extension_rate=0.,
                               extension_err=5.,
                               dextral_rate=0.,
                               dextral_err=5.,
                               lsd=20.,
                               usd=0.,
                               name="",
                               hw="1",
                               fw="3",
                               fid=4
                               )

    faults_ans = convert(Vector{Oiler.Fault}, [fault_ans_1; fault_ans_2; 
                                               fault_ans_3])

    for (i, fault) in enumerate(faults)
        @test struct_equals(fault, faults_ans[i])
    end

end


function test_process_faults_from_gis_files_onefile()
    fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
        test_gpkg; layernames=["faults"])

    fault_ans_1 = Oiler.Fault(;trace=[0.49262599938633844 0.7071696809078425; 
                               0.6397098710932747 1.7367567828563972],
                               dip=60.,
                               dip_dir="E",
                               extension_rate=0.,
                               extension_err=1.,
                               dextral_rate=0.,
                               dextral_err=1.,
                               lsd=20.,
                               usd=0.,
                               name="f1",
                               hw="3",
                               fw="5",
                               fid=1
                               )

    fault_ans_2 = Oiler.Fault(;trace=[1.6490095424615623 2.233798832072941; 
                                    1.8163808447487657 1.4882357582481254; 
                                    0.9288057568620807 1.2904333100905214],
                               dip=30.,
                               dip_dir="NW",
                               extension_rate=1.,
                               extension_err=5.,
                               dextral_rate=3.,
                               dextral_err=2.,
                               lsd=15.,
                               usd=0.,
                               name="f2",
                               hw="1",
                               fw="2",
                               fid=2
                               )

    fault_ans_3 = Oiler.Fault(;trace=[0.9288057568620807 1.2904333100905214; 
                                    0.6397098710932747 1.7367567828563972],
                               dip=89.,
                               dip_dir="E",
                               extension_rate=0.,
                               extension_err=5.,
                               dextral_rate=0.,
                               dextral_err=5.,
                               lsd=20.,
                               usd=0.,
                               name="",
                               hw="1",
                               fw="3",
                               fid=4
                               )

    faults_ans = convert(Vector{Oiler.Fault}, [fault_ans_1; fault_ans_2; 
                                               fault_ans_3])

    for (i, fault) in enumerate(faults)
        @test struct_equals(fault, faults_ans[i])
    end

    for (i, vel) in enumerate(fault_vels)
    end

end


function test_make_vels_from_faults()
    fault_df = Oiler.IO.gis_vec_file_to_df(test_gpkg; layername="faults")
    faults = Oiler.IO.process_faults_from_df(fault_df)
    fault_vels = Oiler.IO.make_vels_from_faults(faults)

    fault_vels_ans = [
        Oiler.VelocityVectorSphere(0.5171328813287053, 0.8787680378934949, 
        0.0, 0.0, 0.0, 0.9999999999999999, 0.9999999999999999, 0.0, 0.0, 
        "3", "5", "f1", "fault");
        Oiler.VelocityVectorSphere(0.5661538404433727, 1.2219642381684521, 
        0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, "3", "5", "f1", "fault");
        Oiler.VelocityVectorSphere(0.6151873278949507, 1.5651595439933637, 0.0, 
        0.0, 0.0, 1.0, 1.0, 0.0, 0.0, "3", "5", "f1", "fault");
        Oiler.VelocityVectorSphere(1.5289218882820936, 2.076591720553694, 
        -1.025534594670751, -2.991367378830875, 0.0, 4.155447322858261, 
        3.4252383197304845, 0.0, -10.128652512234064, "1", "2", "f2", "fault");
        Oiler.VelocityVectorSphere(1.288816444978112, 1.7621508521794715, 
        -1.025080257472689, -2.9915231013214196, 0.0, 4.155817486495381, 
        3.4247891933576313, 0.0, -10.127811322007991, "1", "2", "f2", "fault");
        Oiler.VelocityVectorSphere(1.0487920710603735, 1.4476790679944265, 
        -1.0246948163748018, -2.9916551494606813, 0.0, 4.156131453124761, 
        3.424408174319625, 0.0, -10.127096992327228, "1", "2", "f2", "fault");
        Oiler.VelocityVectorSphere(0.8565426682044769, 1.4020177003206074, 0.0, 
        0.0, 0.0, 5.0, 5.0, 0.0, 0.0, "1", "3", "", "fault");
        Oiler.VelocityVectorSphere(0.7119952745185105, 1.6251796144959085, 0.0, 
        0.0, 0.0, 5.0, 5.0, 0.0, 0.0, "1", "3", "", "fault")
    ]

    @test fault_vels == fault_vels_ans
end


function test_make_vel_from_slip_rate()
    slip_rate_df = Oiler.IO.gis_vec_file_to_df(test_gpkg; layername="geol_slip_rates")
    fault_df = Oiler.IO.gis_vec_file_to_df(test_gpkg; layername="faults")

    vel = Oiler.IO.make_vel_from_slip_rate(slip_rate_df[1,:], fault_df)

    vel_ans = Oiler.VelocityVectorSphere(
    0.5858878181236526,
    1.3583165599070586,
    0.3534116980833591,
    2.4748939718011833,
    0.0,
    0.9905882220926491,
    0.28501749814936717,
    0.0,
    -0.1311984721725788,
    "3",
    "5",
    "1",
    "fault"
    )

    @test vel == vel_ans

end


function test_make_geol_slip_rate_vels()
    slip_rate_df = Oiler.IO.gis_vec_file_to_df(test_gpkg; layername="geol_slip_rates")
    fault_df = Oiler.IO.gis_vec_file_to_df(test_gpkg; layername="faults")

    vels = Oiler.IO.make_geol_slip_rate_vels(slip_rate_df, fault_df)
    vels_ans = [VelocityVectorSphere(0.5858878181236526, 1.3583165599070586, 
        0.3534116980833591, 2.4748939718011833, 0.0, 0.9905882220926491, 
        0.28501749814936717, 0.0, -0.1311984721725788, "3", "5", "1", "fault");
        VelocityVectorSphere(0.8377540074757548, 1.4320519252424566, 
        -0.1678871580639106, -0.10869177594106878, 0.0, 0.54990369919767, 
        0.8411931535674303, 0.0, -0.45163684584971786, "1", "3", "4", "fault");
    VelocityVectorSphere(1.482047800040411, 1.4147786332165946, 
    0.07204703826822562, -0.21168189406932478, 0.0, 0.01587110108066162, 
    0.015751449155154633, 0.0, -0.00014998806921026007, "1", "2", "2", "fault")]

    @test vels == vels_ans

end


function test_get_coords_from_geom_polyline()
    fault_df, block_df, gnss_df = load_geodataframes()
    trace = Oiler.IO.get_coords_from_geom(fault_df[2, :geometry])
    trace_ans = [1.64901   2.2338;
                 1.81638   1.48824;
                 0.928806  1.29043]
    @test isapprox(trace, trace_ans, rtol=0.001)
end


function test_get_block_idx_for_points()
    fault_df, block_df, gnss_df = load_geodataframes()
    bx = Oiler.IO.get_block_idx_for_points(gnss_df, block_df)
    @test ismissing(bx[1])
    @test bx[2] == 5
    @test bx[3] == 5
    @test bx[4] == 1
end

@testset "io.jl unit tests" begin
    test_get_coords_from_geom_polyline()
    test_gis_vec_file_to_df_gpkg_faults()
    test_get_block_idx_for_points()
    test_gis_vec_file_to_df_geojson_blocks()
    # load_geodataframes
    test_row_to_fault()
    test_process_faults_from_df()
    test_process_faults_from_gis_files_onefile()
    test_make_vels_from_faults()
    test_make_vel_from_slip_rate()
    test_make_geol_slip_rate_vels()
end
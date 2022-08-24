using Test

using Oiler


# from https://stackoverflow.com/questions/66801702/
# deriving-equality-for-julia-structs-with-mutable-members
function struct_equals(a::S, b::S) where {S}
    for name in fieldnames(S)
        val1 = getfield(a, name)
        val2 = getfield(b, name)
        try
            if !isapprox(val1, val2; atol=1e-12)
                println(name, ": ", val1, " != ", val2)
                return false
            end
        catch MethodError
            if val1 != val2
                println(name, ": ", val1, " != ", val2)
                return false
            end
        end
    end
    return true
end


test_gpkg = "./test_data/io_test.gpkg"
test_gj_blocks = "./test_data/io_test_blocks.geojson"
test_gj_faults = "./test_data/io_test_faults.geojson"
test_gj_boundary = "./test_data/io_test_boundary.geojson"
test_gj_gnss_vels = "./test_data/io_test_gnss_vels.geojson"
test_gj_geol_slip_rates = "./test_data/io_test_geol_slip_rates.geojson"


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
    block_df = Oiler.IO.gis_vec_file_to_df(test_gj_blocks)
    @test size(block_df, 1) == 4
end


function load_geodataframes()
    fault_df = Oiler.IO.gis_vec_file_to_df(test_gj_faults)
    block_df = Oiler.IO.gis_vec_file_to_df(test_gj_blocks)
    gnss_df = Oiler.IO.gis_vec_file_to_df(test_gj_gnss_vels)

    (fault_df, block_df, gnss_df)
end


function test_row_to_fault()
    fault_df = Oiler.IO.gis_vec_file_to_df(test_gj_faults)

    fault = Oiler.IO.row_to_fault(fault_df[1, :])

    fault_ans = Oiler.Fault(; trace=[0.492625999386338 0.707169680907842
            0.639709871093275 1.736756782856397],
        dip=60.0,
        dip_dir="E",
        extension_rate=0.0,
        extension_err=1.0,
        dextral_rate=0.0,
        dextral_err=1.0,
        lsd=20.0,
        usd=0.0,
        name="f1",
        hw="3",
        fw="5",
        fid=1
    )

    @test struct_equals(fault, fault_ans)
end


function test_process_faults_from_df()

    fault_df = Oiler.IO.gis_vec_file_to_df(test_gj_faults)

    faults = Oiler.IO.process_faults_from_df(fault_df)

    fault_ans_1 = Oiler.Fault(; trace=[0.492625999386338 0.707169680907842
            0.639709871093275 1.736756782856397],
        dip=60.0,
        dip_dir="E",
        extension_rate=0.0,
        extension_err=1.0,
        dextral_rate=0.0,
        dextral_err=1.0,
        lsd=20.0,
        usd=0.0,
        name="f1",
        hw="3",
        fw="5",
        fid=1
    )

    fault_ans_2 = Oiler.Fault(; trace=[1.649009542461562 2.233798832072941
            1.816380844748766 1.488235758248125
            0.928805756862081 1.290433310090521],
        dip=30.0,
        dip_dir="NW",
        extension_rate=1.0,
        extension_err=5.0,
        dextral_rate=3.0,
        dextral_err=2.0,
        lsd=15.0,
        usd=0.0,
        name="f2",
        hw="1",
        fw="2",
        fid=2
    )

    fault_ans_3 = Oiler.Fault(; trace=[0.928805756862081 1.290433310090521
            0.639709871093275 1.736756782856397],
        dip=89.0,
        dip_dir="E",
        extension_rate=0.0,
        extension_err=5.0,
        dextral_rate=0.0,
        dextral_err=5.0,
        lsd=20.0,
        usd=0.0,
        name="",
        hw="1",
        fw="3",
        fid=4
    )

    faults_ans = [fault_ans_1 fault_ans_2 fault_ans_3]

    for (i, fault) in enumerate(faults)
        @test struct_equals(fault, faults_ans[i])
    end

end


function test_process_faults_from_gis_files_onefile()
    fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
        test_gj_faults)

    fault_ans_1 = Oiler.Fault(; trace=[0.492625999386338 0.707169680907842
            0.639709871093275 1.736756782856397],
        dip=60.0,
        dip_dir="E",
        extension_rate=0.0,
        extension_err=1.0,
        dextral_rate=0.0,
        dextral_err=1.0,
        lsd=20.0,
        usd=0.0,
        name="f1",
        hw="3",
        fw="5",
        fid=1
    )

    fault_ans_2 = Oiler.Fault(; trace=[1.649009542461562 2.233798832072941
            1.816380844748766 1.488235758248125
            0.928805756862081 1.290433310090521],
        dip=30.0,
        dip_dir="NW",
        extension_rate=1.0,
        extension_err=5.0,
        dextral_rate=3.0,
        dextral_err=2.0,
        lsd=15.0,
        usd=0.0,
        name="f2",
        hw="1",
        fw="2",
        fid=2
    )

    fault_ans_3 = Oiler.Fault(; trace=[0.928805756862081 1.290433310090521
            0.639709871093275 1.736756782856397],
        dip=89.0,
        dip_dir="E",
        extension_rate=0.0,
        extension_err=5.0,
        dextral_rate=0.0,
        dextral_err=5.0,
        lsd=20.0,
        usd=0.0,
        name="",
        hw="1",
        fw="3",
        fid=4
    )

    faults_ans = [fault_ans_1 fault_ans_2 fault_ans_3]

    for (i, fault) in enumerate(faults)
        @test struct_equals(fault, faults_ans[i])
    end

end


function test_get_blocks_in_bounds()
    block_df = Oiler.IO.gis_vec_file_to_df(test_gj_blocks)
    bound_df = Oiler.IO.gis_vec_file_to_df(test_gj_boundary)

    block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df)

    @test block_df[:, :fid] == [1, 3, 5]
end


function test_process_faults_from_gis_files_w_bounds()

    block_df = Oiler.IO.gis_vec_file_to_df(test_gj_blocks)
    bound_df = Oiler.IO.gis_vec_file_to_df(test_gj_boundary)
    block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df)

    fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
        test_gj_faults, block_df=block_df,
        subset_in_bounds=true)

    fault_ans_1 = Oiler.Fault(; trace=[0.492625999386338 0.707169680907842
            0.639709871093275 1.736756782856397],
        dip=60.0,
        dip_dir="E",
        extension_rate=0.0,
        extension_err=1.0,
        dextral_rate=0.0,
        dextral_err=1.0,
        lsd=20.0,
        usd=0.0,
        name="f1",
        hw="3",
        fw="5",
        fid=1
    )

    fault_ans_3 = Oiler.Fault(; trace=[0.928805756862081 1.290433310090521
            0.639709871093275 1.736756782856397],
        dip=89.0,
        dip_dir="E",
        extension_rate=0.0,
        extension_err=5.0,
        dextral_rate=0.0,
        dextral_err=5.0,
        lsd=20.0,
        usd=0.0,
        name="",
        hw="1",
        fw="3",
        fid=4
    )

    @test struct_equals(faults[1], fault_ans_1)
    @test struct_equals(faults[2], fault_ans_3)
end

function test_make_vels_from_faults()
    fault_df = Oiler.IO.gis_vec_file_to_df(test_gj_faults)
    faults = Oiler.IO.process_faults_from_df(fault_df)
    fault_vels = Oiler.IO.make_vels_from_faults(faults)

    fault_vels_ans = [
        Oiler.VelocityVectorSphere(0.5073298829822619, 0.810128712631411,
            0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, "3", "5", "f1", "fault"),
        Oiler.VelocityVectorSphere(0.5367399870040994, 1.0160466092429714,
            0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, "3", "5", "f1", "fault"),
        Oiler.VelocityVectorSphere(0.5661538404433727, 1.2219642381684515,
            0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, "3", "5", "f1", "fault"),
        Oiler.VelocityVectorSphere(0.5955722040155479, 1.4278815450968356,
            0.0, 0.0, 0.0, 0.9999999999999999, 0.9999999999999999,
            0.0, 0.0, "3", "5", "f1", "fault"),
        Oiler.VelocityVectorSphere(0.6249958388976995, 1.6337984756838961,
            0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, "3", "5", "f1", "fault"),
        Oiler.VelocityVectorSphere(1.5889626154035639, 2.15519645864291,
            -1.0255345946707481, -2.9913673788308754, 0.0, 4.155447322858263,
            3.425238319730482, 0.0, -10.128652512234057, "1", "2", "f2", "fault"),
        Oiler.VelocityVectorSphere(1.4688871343408456, 1.997984704357551,
            -1.0252988111264187, -2.991448202443551, 0.0, 4.155639433917717,
            3.425005240181106, 0.0, -10.128216078017015,
            "1", "2", "f2", "fault"),
        Oiler.VelocityVectorSphere(1.3488346392904595, 1.8407641836954147,
            -1.0250802574726836, -2.991523101321422, 0.0, 4.155817486495385,
            3.424789193357626, 0.0, -10.127811322007982, "1", "2", "f2", "fault"),
        Oiler.VelocityVectorSphere(1.2288033174347306, 1.6835355884591243,
            -1.024878927726256, -2.991592081735523, 0.0, 4.155981490007017,
            3.4245901732497948, 0.0, -10.127438281884793, "1", "2", "f2", "fault"),
        Oiler.VelocityVectorSphere(1.1087913568609091, 1.5262996099996498,
            -1.0246948163747942, -2.991655149460684, 0.0, 4.156131453124767,
            3.424408174319617, 0.0, -10.12709699232721, "1", "2", "f2", "fault"),
        Oiler.VelocityVectorSphere(0.9887969464837787, 1.369056939255023,
            -1.0245279183769676, -2.9917123097761515, 0.0, 4.15626738377776,
            3.4242431915016156, 0.0, -10.126787485025112, "1", "2", "f2", "fault"),
        Oiler.VelocityVectorSphere(0.8806311027661972, 1.3648231425242108,
            0.0, 0.0, 0.0, 5.0, 5.0, 0.0, 0.0, "1", "3", "", "fault"),
        Oiler.VelocityVectorSphere(0.7842726904596813, 1.513599861041644,
            0.0, 0.0, 0.0, 5.0, 5.0, 0.0, -1.7763568394002505e-15,
            "1", "3", "", "fault"),
        Oiler.VelocityVectorSphere(0.6879010546723754, 1.6623722999657207,
            0.0, 0.0, 0.0, 5.0, 5.0, 0.0, 1.7763568394002505e-15,
            "1", "3", "", "fault")]

    for (i, fault) in enumerate(fault_vels)
        @test struct_equals(fault, fault_vels_ans[i])
    end
end


function test_make_vel_from_slip_rate()
    slip_rate_df = Oiler.IO.gis_vec_file_to_df(test_gj_geol_slip_rates)
    fault_df = Oiler.IO.gis_vec_file_to_df(test_gj_faults)

    vel = Oiler.IO.make_vel_from_slip_rate(slip_rate_df[1, :], fault_df)

    vel_ans = Oiler.VelocityVectorSphere(
        0.585887818123653,
        1.358316559907059,
        0.3534116980833613,
        2.4748939718011833,
        0.0,
        0.990588222092649,
        0.28501749814936755,
        0.0,
        -0.1311984721725796,
        "3",
        "5",
        "1",
        "fault"
    )

    @test vel == vel_ans

end


function test_make_geol_slip_rate_vels()
    slip_rate_df = Oiler.IO.gis_vec_file_to_df(test_gj_geol_slip_rates)
    fault_df = Oiler.IO.gis_vec_file_to_df(test_gj_faults)

    vels = Oiler.IO.make_geol_slip_rate_vels(slip_rate_df, fault_df)
    vels_ans = [VelocityVectorSphere(0.837754007475755, 1.432051925242457,
            -0.1678871580639106, -0.10869177594106878, 0.0, 0.54990369919767,
            0.8411931535674303, 0.0, -0.45163684584971814, "1", "3", "4", "fault"),
        VelocityVectorSphere(1.482047800040411, 1.414778633216595,
            0.09829170611540894, -0.20084506593123488, 0.0, 0.017015853616702464,
            0.014507264583475196, 0.0, -0.00014469500954570033, "1", "2", "2", "fault")]

    for (i, vel) in enumerate(vels)
        @test struct_equals(vel, vels_ans[i])
    end

end


function test_get_coords_from_geom_polyline()
    fault_df, block_df, gnss_df = load_geodataframes()
    trace = Oiler.IO.get_coords_from_geom(fault_df[2, :geometry])
    trace_ans = [1.64901 2.2338
        1.81638 1.48824
        0.928806 1.29043]
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
    ##test_gis_vec_file_to_df_gpkg_faults()
    test_get_block_idx_for_points()
    test_gis_vec_file_to_df_geojson_blocks()
    ## load_geodataframes
    test_row_to_fault()
    test_process_faults_from_df()
    test_get_blocks_in_bounds()
    test_process_faults_from_gis_files_onefile()
    test_process_faults_from_gis_files_w_bounds()
    test_make_vels_from_faults()
    test_make_vel_from_slip_rate()
    test_make_geol_slip_rate_vels()
end

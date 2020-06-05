using Revise

using Test

using Oiler


test_gpkg = "./test_data/io_test.gpkg"

function test_gis_vec_file_to_df_gpkg_faults()
    fault_df = Oiler.IO.gis_vec_file_to_df(test_gpkg;
        layername="faults")

    @test size(fault_df) == (3,15)
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

fault_df = Oiler.IO.gis_vec_file_to_df(test_gpkg; layername="faults")
block_df = Oiler.IO.gis_vec_file_to_df(test_gpkg; layername="blocks")
gnss_df = Oiler.IO.gis_vec_file_to_df(test_gpkg;layername="gnss_vels")


function test_get_coords_from_geom_polyline()
    trace = Oiler.IO.get_coords_from_geom(fault_df[2, :geometry])
    trace_ans = [1.64901   2.2338;
                 1.81638   1.48824;
                 0.928806  1.29043]
    @test isapprox(trace, trace_ans, rtol=0.001)
end


function test_get_block_idx_for_points()
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
end
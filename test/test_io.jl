using Revise

using Oiler

test_gpkg = "./test_data/io_test.gpkg"

function test_gis_vec_file_to_df_gpkg()
    line_df = Oiler.IO.gis_vec_file_to_df(test_gpkg;
        layername="line_test")
end

line_df = test_gis_vec_file_to_df_gpkg()
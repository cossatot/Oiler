using Test

using CSV
using JSON
# using PyPlot
using DataFrames

using Oiler

atol = 1e-10

# load tris
tri_json = JSON.parsefile("./test_data/fake_na_ca/ca_na_fake_la_tris.geojson",
    dicttype=Dict)

function tri_from_geojson(tri_feature)
    coords = tri_feature["geometry"]["coordinates"][1]
    Oiler.Tris.Tri(p1=Float64.(coords[1]), 
                   p2=Float64.(coords[2]), 
                   p3=Float64.(coords[3]), 
                   name=tri_feature["properties"]["fid"])
end

tris = map(tri_from_geojson, tri_json["features"])

# load vel points
vels = CSV.read("./test_data/fake_na_ca/ca_na_fake_pts.csv")
lons = [vels[i,:].X for i in 1:size(vels, 1)]
lats = [vels[i,:].Y for i in 1:size(vels, 1)]



function calc_tri_disps(tris, gnss_lons, gnss_lats)

    tri_gnss_partials = hcat( collect(
        [Oiler.Elastic.arrange_tri_partials(
            Oiler.Elastic.calc_tri_slip(tri, gnss_lons, gnss_lats; ds_slip=100.)...)
         for tri in tris]
    )...)

end


tri_ds = calc_tri_disps(tris, lons, lats)
sum_tri_ds = sum(tri_ds; dims=2)

tot_ve = sum_tri_ds[1:3:end]
tot_vn = sum_tri_ds[2:3:end]
tot_vu = sum_tri_ds[3:3:end]

tri_partials = Oiler.Elastic.calc_tri_effects(tris, lons, lats)

tri_slip_inv = tri_partials \ sum_tri_ds

function test_results()
    for i in length(tris)
        @test isapprox(tri_slip_inv[i * 2 - 1], 100.; atol=atol)
        @test isapprox(tri_slip_inv[i * 2], 0.; atol=atol)
    end
end

@testset "test solve only tris" begin
    test_results()
end
using DataFrames
using CSV
using PyPlot


vels = CSV.read("./data/fault_vels.csv");

figure()
quiver(vels.lon, vels.lat, vels.ve, vels.vn, scale = 100)
show()
using CSV
using DataFrames

include("./block_rotations.jl")

function vel_from_row(row::DataFrameRow)
    VelocityVectorSphere(lond = row.lon, latd = row.lat, 
                         ve = row.ve, vn = row.vn, fix = row.fix, mov = row.mov)
end


function load_vels_from_csv(filepath)
    # will need, eventually, to support velocity weights
    vels = CSV.read(filepath)

    vel_array = [vel_from_row(vels[i,:]) for i in 1:size(vels, 1)]
end


function group_vels_by_fix_mov(vels::Array{VelocityVectorSphere})
    # flipping opposite vels for now; maybe reconsider later

    vel_groups = Dict{Tuple{String,String},Array{VelocityVectorSphere,1}}()

    for vel in vels
        fm = (vel.fix, vel.mov)
        if haskey(vel_groups, fm)
            push!(vel_groups[fm], vel)
        elseif haskey(vel_groups, Base.reverse(fm))
            push!(vel_groups[Base.reverse(fm)], reverse(vel))
        else
            vel_groups[fm] = [vel]
        end
    end
    vel_groups
end
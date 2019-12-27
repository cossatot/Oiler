module IO

export vel_from_row, load_vels_from_csv, group_vels_by_fix_mov

using ..Oiler: VelocityVectorSphere

using CSV
using DataFrames: DataFrameRow


function vel_from_row(row::DataFrameRow)
    # maybe use dict here to match keys

    if :ee in names(row)
        ee = row.ee
    else
        ee = 0.
    end

    if :en in names(row)
        en = row.en
    else
        en = 0.
    end
    
    VelocityVectorSphere(lond = row.lon, latd = row.lat, ee = ee, en = en,
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

end
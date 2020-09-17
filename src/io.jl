module IO

export vel_from_row, load_vels_from_csv, group_vels_by_fix_mov

using ..Oiler: VelocityVectorSphere, PoleSphere, PoleCart, pole_cart_to_sphere

using Logging

using CSV
using ArchGDAL
using DataFrames: DataFrame, DataFrameRow

const AG = ArchGDAL


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
    
    VelocityVectorSphere(lon=row.lon, lat=row.lat, ee=ee, en=en,
                         ve=row.ve, vn=row.vn, fix=row.fix, mov=row.mov)
end


function load_vels_from_csv(filepath)
    # will need, eventually, to support velocity weights
    vels = CSV.read(filepath)

    vel_array = [vel_from_row(vels[i,:]) for i in 1:size(vels, 1)]
end

"""
    group_vels_by_fix_mov

Groups velocities by the poles which they describe the relative motion.

# Arguments
- `vels`: An array of `VelocityVectorSphere`s

# Returns
- `vel_groups`: A dictionary with keys representing the poles as
  tuples of format (`fix`, `mov`), and values of arrays of 
  `VelocityVectorSphere` of each of the velocities that shares this pole
  configuration.
"""
function group_vels_by_fix_mov(vels::Array{VelocityVectorSphere})
    # flipping opposite vels for now; maybe reconsider later

    vel_groups = Dict{Tuple{String,String},Array{VelocityVectorSphere,1}}()

    for vel in vels
        fm = (vel.fix, vel.mov)
        if vel.fix == vel.mov
            if vel.vel_type == "GNSS"
                if haskey(vel_groups, fm)
                    push!(vel_groups[fm], vel)
                else
                    vel_groups[fm] = [vel]
                end
            else
                @warn "$vel has same fix, mov, leaving it out."
            end
        elseif haskey(vel_groups, fm)
            push!(vel_groups[fm], vel)
        elseif haskey(vel_groups, Base.reverse(fm))
            push!(vel_groups[Base.reverse(fm)], reverse(vel))
        else
            vel_groups[fm] = [vel]
        end
    end
    vel_groups
end


function pole_to_dict(pole::PoleCart)
    Dict("x" => pole.x, "y" => pole.y, "z" => pole.z, "fix" => pole.fix, 
         "mov" => pole.mov)
end


function pole_to_dict(pole::PoleSphere)
    Dict("lon" => pole.lon, "lat" => pole.lat, "rotrate" => pole.rotrate, 
        "fix" => pole.fix, "mov" => pole.mov)
end


# function pole_dict_to_df(poles::Dict)
#

# end

function poles_to_df(poles::Array{PoleSphere,1})
    df = DataFrame()
    df.lon = [p.lon for p in poles]
    df.lat = [p.lat for p in poles]
    df.rotrate = [p.rotrate for p in poles]
    df.fix = [p.fix for p in poles]
    df.mov = [p.mov for p in poles]

    df
end


function poles_to_df(poles::Array{PoleCart,1};
                     convert_to_sphere::Bool=false)

    if convert_to_sphere
        pole_sphere = map(pole_cart_to_sphere, poles)
        df = poles_to_df(pole_sphere)
    else
        df = DataFrame()
        df.x = [p.x for p in poles]
        df.y = [p.y for p in poles]
        df.z = [p.z for p in poles]
        df.fix = [p.fix for p in poles]
        df.mov = [p.mov for p in poles]
    end
    df
end



function gis_vec_file_to_df(filename::AbstractString; layername="")
    dataset = AG.read(filename)
    if layername == ""
        layer = AG.getlayer(dataset, 0)
    else
        layer = AG.getlayer(dataset, layername)
    end

    nfeat = AG.nfeature(layer)
    nfield = AG.nfield(layer)

    # prepare Dict with empty vectors of the right type for each field
    d = Dict{String,Vector}()
    featuredefn = AG.layerdefn(layer)
    for field_no in 0:nfield - 1
        field = AG.getfielddefn(featuredefn, field_no)
        name = AG.getname(field)
        typ = AG._FIELDTYPE[AG.gettype(field)]
        d[name] = typ[]
    end
    d["geometry"] = AG.IGeometry[]
    d["fid"] = Int64[]

    # loop over the features to fill the vectors in the Dict
    # for fid in 0:nfeat-1
    for nf in 1:nfeat
        try
            feature = AG.unsafe_nextfeature(layer)# do feature
            for (k, v) in pairs(d)
                if k == "geometry"
                    val = AG.getgeom(feature, 0)
                elseif k == "fid"
                    try
                        val = AG.getfield(feature, k)
                    catch
                        val = AG.getfid(feature)
                    end
                else
                    val = AG.getfield(feature, k)
                end
                push!(v, val)
            end
        catch e
            println(e)
        end
    end
    # construct a DataFrame from the Dict
    df = DataFrame(d)
end


function get_geom_coords_from_feature(feature)
    geom = AG.getgeom(feature, 0)
    coords = get_coords_from_geom(geom)
end


function get_coords_from_geom(geom)
    n_pts = AG.ngeom(geom)

    coords = Array{Float64}(undef, n_pts, 2)

    for i in 0:n_pts - 1
        coords[i + 1, 1] = AG.getx(geom, i)
        coords[i + 1, 2] = AG.gety(geom, i)
    end
    coords
end


function get_block_idx_for_points(point_df, block_df)
    point_geoms = point_df[:, :geometry]
    idxs = Array{Any}(missing, length(point_geoms))

    for i_b in 1:size(block_df, 1)
        block_geom = block_df[i_b, :geometry]
        block_fid = block_df[i_b, :fid]
        
        for (i_p, pg) in enumerate(point_geoms)
            if AG.contains(block_geom, pg)
                if !ismissing(idxs[i_p])
                    pfid = point_df[i_p, :fid]
                    @warn "$pfid in multiple blocks"
                else
                    idxs[i_p] = block_fid
                end
            end
        end
    end
    idxs
end





end

module IO

export vel_from_row, load_vels_from_csv, group_vels_by_fix_mov

using ..Oiler
using ..Oiler: VelocityVectorSphere, PoleSphere, PoleCart, pole_cart_to_sphere


using Logging
import Base.Threads.@threads

using CSV
using ArchGDAL
using GeoFormatTypes
using DataFrames: DataFrame, DataFrameRow
using DataFramesMeta

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
    # d["fid"] = Int64[]
    d["fid"] = []

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

    @threads for i_b in 1:size(block_df, 1)
        block_geom = block_df[i_b, :geometry]
        block_fid = block_df[i_b, :fid]
        
        @threads for (i_p, pg) in collect(enumerate(point_geoms))
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


function get_block_idx_for_points(point_df, block_df, to_epsg)

    point_geoms = point_df[:, :geometry]
    point_geoms_p = [AG.clone(p) for p in point_geoms]
    point_geoms_p = [AG.reproject(pg, 
            ProjString("+proj=longlat +datum=WGS84 +no_defs"), 
            EPSG(to_epsg)) for pg in point_geoms_p]

    idxs = Array{Any}(missing, length(point_geoms))

    @threads for i_b in 1:size(block_df, 1)
        block_geom = block_df[i_b, :geometry]

        block_geom_p = AG.reproject(AG.clone(block_geom), 
            ProjString("+proj=longlat +datum=WGS84 +no_defs"), 
            EPSG(to_epsg))

        block_fid = block_df[i_b, :fid]
        
        @threads for (i_p, pg) in collect(enumerate(point_geoms_p))
            if AG.contains(block_geom_p, pg)
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





function val_nothing_fix(vel; return_val=0.)
    if vel == ""
        return return_val
    elseif isnothing(vel) | ismissing(vel)
        return return_val
    else
        if typeof(vel) == String
            return parse(Float64, vel)
        else
            return vel
        end
    end
end


function err_nothing_fix(err; return_val=1., weight=1.)
    if err == ""
        return return_val
    elseif isnothing(err) | ismissing(err) | (err == 0.)
        return return_val
    else
        if typeof(err) == String
            err = parse(Float64, err) / weight
            if iszero(err)
                return return_val
            else
                return err
            end
        else
            return err
        end
    end
end


function row_to_fault(row; name="", dip_dir=:dip_dir, v_ex=:v_ex, e_ex=:e_ex,
        v_rl=:v_rl, e_rl=:e_rl, dip=:dip, hw=:hw, fw=:fw, usd=:usd, lsd=:lsd,
        v_default=0., e_default=5., usd_default=0., lsd_default=20.)

    trace = Oiler.IO.get_coords_from_geom(row[:geometry])
    if name in names(row)
        _name = row[:name]
    else
        _name = name
    end

    Oiler.Fault(trace=trace,
        dip_dir=row[dip_dir],
        extension_rate=val_nothing_fix(row[v_ex], return_val=v_default),
        extension_err=err_nothing_fix(row[e_ex], return_val=e_default),
        dextral_rate=val_nothing_fix(row[v_rl], return_val=v_default),
        dextral_err=err_nothing_fix(row[e_rl], return_val=e_default),
        dip=row[dip],
        name=_name,
        hw=row[hw],
        fw=row[fw],
        usd=val_nothing_fix(row[usd], return_val=usd_default),
        lsd=val_nothing_fix(row[lsd], return_val=lsd_default),
        )
end


function make_vel_from_slip_rate(slip_rate_row, fault_df; err_return_val=1., weight=1.)
    fault_seg = slip_rate_row[:fault_seg]
    fault_idx = fault_seg# parse(Int, fault_seg)
    fault_row = @where(fault_df, :fid .== fault_idx)[1,:]
    fault = row_to_fault(fault_row)

    extension_rate = val_nothing_fix(slip_rate_row[:extension_rate])
    extension_err = err_nothing_fix(slip_rate_row[:extension_err]; 
        return_val=err_return_val, weight=weight)
    dextral_rate = val_nothing_fix(slip_rate_row[:dextral_rate])
    dextral_err = err_nothing_fix(slip_rate_row[:dextral_err]; 
        return_val=err_return_val, weight=weight)

    ve, vn = Oiler.Faults.fault_slip_rate_to_ve_vn(dextral_rate, 
                                                   extension_rate,
                                                   fault.strike)

    ee, en, cen = Oiler.Faults.fault_slip_rate_err_to_ee_en(dextral_err, 
                                                            extension_err,
                                                            fault.strike)

    pt = Oiler.IO.get_coords_from_geom(slip_rate_row[:geometry])
    lon = pt[1]
    lat = pt[2]
    
    VelocityVectorSphere(lon=lon, lat=lat, ve=ve, vn=vn, fix=fault.hw,
                         ee=ee, en=en, cen=cen,
                         mov=fault.fw, vel_type="fault", name=fault_seg) 
end


function make_geol_slip_rate_vel_vec(geol_slip_rate_df, fault_df;
        err_return_val=1., weight=1.)
    geol_slip_rate_vels = []
    for i in 1:size(geol_slip_rate_df, 1)
        slip_rate_row = geol_slip_rate_df[i,:]
        if (slip_rate_row[:include] == true) | (slip_rate_row[:include] == "1")
            push!(geol_slip_rate_vels, make_vel_from_slip_rate(slip_rate_row, 
                                                               fault_df;
                                                               err_return_val=err_return_val,
                                                               weight=weight))
        end
    end

    geol_slip_rate_vels = convert(Array{VelocityVectorSphere}, geol_slip_rate_vels)
end


function tri_from_feature(feat)
    # TODO: add support for reading in rates
    tri = Oiler.Tris.Tri(
       p1=Float64.(feat["geometry"]["coordinates"][1][1]),
       p2=Float64.(feat["geometry"]["coordinates"][1][2]),
       p3=Float64.(feat["geometry"]["coordinates"][1][3]),
       name=string(feat["properties"]["fid"])
       )
end


function tris_from_geojson(tri_json)
    tris = map(tri_from_feature, tri_json["features"])
end


function tri_to_feature(tri)

    gj_feature = Dict(
        "type" => "Feature",
        "geometry" => Dict(
            "type":"Polygon",
            "coordinates" => [[p1, p2, p3, p1]]
        ),
        "properties" => Dict(
            "dip_slip_rate" => round.(tri.dip_slip_rate, digits=3),
            "dip_slip_err" => round.(tri.dip_slip_err, digits=3),
            "strike_slip_rate" => round.(tri.strike_slip_rate, digits=3),
            "strike_slip_err" => round.(tri.strike_slip_err, digits=3),
            "cde" => round.(tri.cde, digits=3),
            "fid" => tri.name,
        )
    )

end

function tris_to_geojson(tris; name="")
    gj = Dict(
        "type" => "FeatureCollection",
        "name":name,
        "features":[tri_to_feature(tri) for tri in tris],
    )
end


function gnss_vel_from_row(row, block; ve=:ve, vn=:vn, ee=:ee, en=:en,
                           fix="1111", name=:station)
    pt = Oiler.IO.get_coords_from_geom(row[:geometry])
    lon = pt[1]
    lat = pt[2]
    Oiler.VelocityVectorSphere(lon=lon, lat=lat, ve=row[ve],
        vn=row[vn], ee=row[ee], en=row[en], name=row[name],
        fix=fix, mov=string(block), vel_type="GNSS")
end


function make_vels_from_gnss_and_blocks(gnss_df, block_df; ve=:ve, vn=:vn, ee=:ee, en=:en,
    fix="1111", name=:station, epsg=0)

    if epsg != 0
        block_idx = get_block_idx_for_points(gnss_df, block_df, epsg)
    else
        block_idx = get_block_idx_for_points(gnss_df, block_df)
    end

    gnss_vels = []

    for (i, block) in enumerate(block_idx)
        if !ismissing(block)
            gv = gnss_vel_from_row(gnss_df[i,:], block; ve=ve, vn=vn, ee=ee,
                    en=en, fix=fix, name=name)
            push!(gnss_vels, gv)
        end
    end

    gnss_vels = convert(Array{VelocityVectorSphere}, gnss_vels)
end



function write_faults_with_rates(faults, filename)
end




end

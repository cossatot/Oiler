module IO

export vel_from_row, load_vels_from_csv, group_vels_by_fix_mov

using ..Oiler
using ..Oiler: VelocityVectorSphere, PoleSphere, PoleCart, pole_cart_to_sphere


using Logging
import Base.Threads.@threads

using CSV
using JSON
using Proj: Transformation
using PolygonOps: inpolygon
using DataFrames: DataFrame, DataFrameRow
using DataFramesMeta



function vel_from_row(row::DataFrameRow)
    # maybe use dict here to match keys

    if :ve in names(row)
        ve = row.ve
    else
        ve = 0.0
    end

    if :vn in names(row)
        vn = row.vn
    else
        vn = 0.0
    end

    if :ee in names(row)
        ee = row.ee
    else
        ee = 0.0
    end

    if :en in names(row)
        en = row.en
    else
        en = 0.0
    end

    VelocityVectorSphere(lon=row.lon, lat=row.lat, ee=ee, en=en,
        ve=ve, vn=vn, fix=row.fix, mov=row.mov)
end


function load_vels_from_csv(filepath)
    # will need, eventually, to support velocity weights
    vels = CSV.read(filepath, DataFrame)

    vel_array = [vel_from_row(vels[i, :]) for i in 1:size(vels, 1)]
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


function gj_linestring_to_array(gj_coords)
    if length(gj_coords[1]) == 2
        arr = vcat([[c[1] c[2]] for c in gj_coords]...)
    elseif length(gj_coords[1]) == 3
        arr = vcat([[c[1] c[2] c[3]] for c in gj_coords]...)
    end
    arr
end


function gj_polygon_to_array(gj_coords; check_poly_winding_order=true,
    epsg=102016)

    coords = gj_linestring_to_array(gj_coords[1])

    if check_poly_winding_order
        trans = make_trans_from_wgs84(epsg)
        trans_coords = trans.(Oiler.Utils.matrix_to_rows(coords))

        # geometries should be clockwise for block-fault trace ordering
        # clockwise is winding order -1
        if Oiler.Geom.check_winding_order(trans_coords) == 1
            reverse!(coords, dims=1)
        end
    end
    coords
end


function gj_point_to_array(gj_coords)
    c = gj_coords
    if length(gj_coords) == 2
        arr = [c[1] c[2]]
    elseif length(gj_coords) == 3
        arr = [c[1] c[2] c[3]]
    end
    arr
end


function reformat_feature(feat; convert_geom=true, check_poly_winding_order=true)
    ref = Dict{String,Any}("geometry" => feat["geometry"]["coordinates"])
    ref["geom_type"] = feat["geometry"]["type"]
    if convert_geom == true
        if ref["geom_type"] == "LineString"
            ref["geometry"] = Oiler.Geom.LineString(gj_linestring_to_array(ref["geometry"]))
        elseif ref["geom_type"] == "Polygon"
            if haskey(feat["properties"], "fid") && (feat["properties"]["fid"] == "ant")
                ref["geometry"] = Oiler.Geom.Polygon(
                    gj_polygon_to_array(ref["geometry"];
                        check_poly_winding_order=check_poly_winding_order,
                        epsg=102019))
            else
                ref["geometry"] = Oiler.Geom.Polygon(
                    gj_polygon_to_array(ref["geometry"];
                        check_poly_winding_order=check_poly_winding_order))
            end
        elseif ref["geom_type"] == "Point"
            ref["geometry"] = Oiler.Geom.Point(gj_point_to_array(ref["geometry"]))
        end
    end

    for (prop, val) in feat["properties"]
        if !isnothing(val)
            ref[prop] = val
        else
            ref[prop] = missing
        end
    end
    return ref
end


function reformat_geojson(gj; convert_geom=true)
    [reformat_feature(feat; convert_geom=convert_geom) for feat in gj["features"]]
end


function tryget(dict, key; default=missing)
    if haskey(dict, key)
        return dict[key]
    else
        return default
    end
end



function gj_to_df(gj; convert_geom=true, fid_drop=[])
    feats = reformat_geojson(gj; convert_geom=convert_geom)

    cols = collect(keys(feats[1]))
    col_dict = Dict{String,Any}()

    for col in cols
        #col_dict[col] = [f[col] for f in feats]
        col_dict[col] = [tryget(f, col) for f in feats]
    end

    if !("fid" in cols)
        col_dict["fid"] = collect(range(1, length=length(feats)))
    end

    df = DataFrame(col_dict)

    df[!, :fid] = string.(df[!, :fid])

    if fid_drop != []
        df = filter(row -> !(row.fid in fid_drop), df)
    end
    df
end


function gis_vec_file_to_df(filename::AbstractString; fid_drop=[], lon=:lon, lat=:lat)
    if !(last(split(filename, ".")) in ["json", "geojson", "csv"])
        throw(ArgumentError("only geojson and csv files currently implemented"))
    end

    if last(split(filename, ".")) == "csv"
        df_raw = CSV.read(filename, DataFrame)
        return georeference_csv!(df_raw, lon=lon, lat=lat)
    
    else
        gj = JSON.parsefile(filename)
        return gj_to_df(gj; fid_drop=fid_drop)
    end
end


function georeference_csv!(df; lon=:lon, lat=:lat)
    df[!, "geometry"] = [Oiler.Geom.Point([row.lon row.lat]) for row in eachrow(df)]
    df
end


function get_coords_from_geom(geom::Oiler.Geom.Point)
    return geom.coords
end


function get_coords_from_geom(geom::Oiler.Geom.LineString)
    return geom.coords
end

function get_coords_from_geom(geom::Oiler.Geom.Polygon)
    return geom.coords
end



function get_block_idx_for_points(point_df, block_df)
    point_geoms = [p.coords for p in point_df[:, :geometry]]
    idxs = Array{Any}(missing, length(point_geoms))

    @threads for i_b in 1:size(block_df, 1)
        block_geom = collect(eachrow(block_df[i_b, :geometry].coords))
        block_fid = block_df[i_b, :fid]

        @threads for (i_p, pg) in collect(enumerate(point_geoms))
            if inpolygon(pg, block_geom) == 1
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


function make_trans_from_wgs84(to_epsg)
    trans = Transformation("EPSG:4326", "EPSG:4326", always_xy=true)
    try
        trans = Transformation("EPSG:4326", "EPSG:$to_epsg", always_xy=true)
    catch
        trans = Transformation("EPSG:4326", "ESRI:$to_epsg", always_xy=true)
    end
    trans
end


function make_trans_to_wgs84(to_epsg)
    trans = Transformation("EPSG:4326", "EPSG:4326", always_xy=true)
    try
        trans = Transformation("EPSG:$to_epsg", "EPSG:4326", always_xy=true)
    catch
        trans = Transformation("ESRI:$to_epsg", "EPSG:4326", always_xy=true)
    end
    trans
end



function get_block_idx_for_points(point_df, block_df, to_epsg)
    point_geoms = [tuple(p.coords...) for p in point_df[:, :geometry]]

    trans = make_trans_from_wgs84(to_epsg)

    point_geoms = trans.(point_geoms)
    point_geoms = [collect(pg) for pg in point_geoms]

    idxs = Array{Any}(missing, length(point_geoms))

    @threads for i_b in 1:size(block_df, 1)
        block_geom = collect(eachrow(block_df[i_b, :geometry].coords))
        block_geom = trans.(block_geom)
        block_fid = block_df[i_b, :fid]

        @threads for (i_p, pg) in collect(enumerate(point_geoms))
            if inpolygon(pg, block_geom) == 1
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


function check_missing(val)
    if isnothing(val) | ismissing(val) | (typeof(val) == Missing)
        return false
    elseif val == ""
        return false
    elseif val == 0.0
        return false
    else
        return true
    end
end


function val_nothing_fix(vel; return_val=0.0)
    if isnothing(vel) | ismissing(vel) | (typeof(vel) == Missing)
        return return_val
    elseif vel == ""
        return return_val
    else
        if typeof(vel) == String
            return parse(Float64, vel)
        else
            return vel
        end
    end
end


function err_nothing_fix(err; return_val=1.0, weight=1.0)
    if isnothing(err) | ismissing(err) | (err == 0.0)
        return return_val
    elseif err == ""
        return return_val
    else
        if typeof(err) == String
            err = parse(Float64, err) / weight
            if iszero(err)
                return return_val / weight
            else
                return err
            end
        else
            return err
        end
    end
end


function adjust_fault_err_by_dip(dip, ds_err, ss_err, dip_adj_remainder)

    new_ds = (dip_adj_remainder * ds_err) + ((1. - dip_adj_remainder) * ds_err * cosd(dip))
    new_ss = (dip_adj_remainder * ss_err) + ((1. - dip_adj_remainder) * ss_err * sind(dip))

    (new_ds, new_ss)
end

function adjust_fault_err_by_dip_proj(dip, ds_err, ss_err, dip_adj_remainder)


    ds_remain =  dip_adj_remainder * ds_err
    ss_remain =  dip_adj_remainder * ss_err

    ds_mod =  (1. - dip_adj_remainder) * ds_err
    ss_mod =  (1. - dip_adj_remainder) * ss_err

    #ds_proj = ds_mod * cosd(dip) - ss_mod * sind(dip)
    #ss_proj = ds_mod * sind(dip) + ss_mod * cosd(dip)
    
    err_proj = Oiler.Geom.rotate_velocity(ds_mod, ss_mod, deg2rad(dip - 90.))

    ds_proj = err_proj[1]
    ss_proj = err_proj[2]

    new_ds = 0. * ds_remain + ds_proj
    new_ss = 0. * ss_remain + abs(ss_proj)

    (new_ds, new_ss)
end

function row_to_fault(row; name="name", dip_dir=:dip_dir, v_ex=:v_ex, e_ex=:e_ex,
    v_rl=:v_rl, e_rl=:e_rl, dip=:dip, hw=:hw, fw=:fw, usd=:usd, lsd=:lsd,
    v_default=0.0, e_default=5.0, usd_default=0.0, lsd_default=20.0, 
    ridge_usd_default=1.0, ridge_lsd_default=2.0,
    transform_usd_default=1.0, transform_lsd_default=10.0,
    ridge_v_rl=0.0, ridge_e_rl=0.5,
    transform_v_ex=0.0, transform_e_ex=0.5,
    adjust_err_by_dip=false, dip_adj_remainder=0.2, fid=:fid)

    trace = Oiler.IO.get_coords_from_geom(row[:geometry])
    if name in names(row)
        _name = row[:name]
        if ismissing(_name)
            _name = ""
        end
    else
        _name = ""
    end
    
    if adjust_err_by_dip
        dp = row[dip]
        #e_default_ds = (dip_adj_remainder * e_default) + (1-dip_adj_remainder * e_default * cosd(dp))
        #e_default_ss = (dip_adj_remainder * e_default) + (1-dip_adj_remainder * e_default * sind(dp))
        e_default_ds, e_default_ss = adjust_fault_err_by_dip(dp, e_default, e_default, dip_adj_remainder)
    else
        e_default_ds = e_default
        e_default_ss = e_default
    end

    if "type" in names(row)
        if row[:type] == "ridge"
            usd_default = ridge_usd_default
            lsd_default = ridge_lsd_default
            v_rl = ridge_v_rl
            e_rl = ridge_e_rl
        elseif row[:type] == "transform"
            usd_default = transform_usd_default
            lsd_default = transform_lsd_default
            v_ex = transform_v_ex
            v_rl = transform_e_ex
        end
    end



    Oiler.Fault(trace=trace,
        dip_dir=row[dip_dir],
        extension_rate=val_nothing_fix(row[v_ex], return_val=v_default),
        extension_err=err_nothing_fix(row[e_ex], return_val=e_default_ds),
        dextral_rate=val_nothing_fix(row[v_rl], return_val=v_default),
        dextral_err=err_nothing_fix(row[e_rl], return_val=e_default_ss),
        dip=row[dip],
        name=_name,
        hw=row[hw],
        fw=row[fw],
        usd=val_nothing_fix(row[usd], return_val=usd_default),
        lsd=val_nothing_fix(row[lsd], return_val=lsd_default),
        fid=row[fid]
    )

end


function process_faults_from_df(fault_df; name="name", dip_dir=:dip_dir, v_ex=:v_ex,
    e_ex=:e_ex, v_rl=:v_rl, e_rl=:e_rl, dip=:dip,
    hw=:hw, fw=:fw, usd=:usd, lsd=:lsd,
    v_default=0.0, e_default=5.0, usd_default=1.0,
    adjust_err_by_dip=false, dip_adj_remainder=0.1,
    lsd_default=15.0, ridge_usd_default=1.0, ridge_lsd_default=2.0,
    transform_usd_default=1.0, transform_lsd_default=10.0,
    ridge_v_rl=0.0, ridge_e_rl=0.5,
    transform_v_ex=0.0, transform_e_ex=0.5,
    fid=:fid, fid_drop=[])

    faults = []
    #for i in 1:size(fault_df, 1)
    for row in eachrow(fault_df)
        #try
        push!(faults, Oiler.IO.row_to_fault(row;
            name=name,
            dip_dir=dip_dir,
            v_ex=v_ex,
            e_ex=e_ex,
            v_rl=v_rl,
            e_rl=e_rl,
            dip=dip,
            hw=hw,
            fw=fw,
            usd=usd,
            lsd=lsd,
            v_default=v_default,
            e_default=e_default,
            usd_default=usd_default,
            lsd_default=lsd_default,
            ridge_usd_default=ridge_usd_default,
            ridge_lsd_default=ridge_lsd_default,
            transform_usd_default=transform_usd_default,
            transform_lsd_default=transform_lsd_default,
            ridge_v_rl=ridge_v_rl,
            ridge_e_rl=ridge_e_rl,
            transform_v_ex=transform_v_ex,
            transform_e_ex=transform_e_ex,
            adjust_err_by_dip=adjust_err_by_dip,
            dip_adj_remainder=dip_adj_remainder,
            fid=fid))
        #catch
        #    prob_row = fault_df[i,:]
        #    warn_msg = "Problem with $prob_row"
        #    @warn warn_msg 
        #end 
    end

    faults = convert(Array{Oiler.Faults.Fault}, faults)

    faults = filter(x -> x.fw != "", faults)
    faults = filter(x -> x.hw != "", faults)
    faults = filter(x -> !(x.fid in fid_drop), faults)

    faults
end


function make_vels_from_faults(faults)
    fault_vels = reduce(vcat, map(Oiler.fault_to_vels, faults))
    fault_vels
end



function load_faults_from_gis_files(fault_files...; layernames=[])

    if length(fault_files) == 1
        if length(layernames) == 0
            fault_df = gis_vec_file_to_df(fault_files[1])
        else
            fault_df = gis_vec_file_to_df(fault_files[1];
                layername=layernames[1])
        end
    else
        if length(layernames) == 0
            fault_df = vcat([gis_vec_file_to_df(file)
                             for file in fault_files]...;
                cols=:union)
        else
            fault_df = vcat([gis_vec_file_to_df(file; layername=layernames[i])
                             for (i, file) in enumerate(fault_files)]...;
                cols=:union)
        end
    end
    fault_df
end

function get_blocks_in_bounds!(block_df, bound_df; epsg=102016)
    bound_orig = collect(eachrow(bound_df[1, :geometry].coords))
    block_geoms_orig = block_df[:, :geometry]
    block_geoms_orig = [collect(eachrow(geom.coords)) for geom in block_geoms_orig]

    if epsg != 0
        trans = make_trans_from_wgs84(epsg)

        bound = trans.(bound_orig)
        block_geoms = [trans.(g) for g in block_geoms_orig]
    else
        bound = bound_orig
        block_geoms = block_geoms_orig
    end

    block_idxs = falses(length(block_geoms))
    for (i, block_geom) in enumerate(block_geoms)
        if any((inpolygon(pt, bound) == 1) for pt in block_geom) |
           any((inpolygon(pt, block_geom) == 1) for pt in bound)
            block_idxs[i] = true
        end
    end
    block_df[block_idxs, :]
end



function get_blocks_and_faults_in_bounds!(block_df, fault_df, bound_df)
    block_df = get_blocks_in_bounds!(block_df, bound_df)
    fault_df = get_faults_bounding_blocks!(fault_df, block_df)

    block_df, fault_df
end


function get_faults_bounding_blocks!(fault_df, block_df)
    block_fids = string.(block_df[:, :fid])
    fault_idxs = falses(size(fault_df, 1))
    for i in 1:length(fault_idxs)
        if (fault_df[i, :hw] in block_fids) & (fault_df[i, :fw] in block_fids)
            fault_idxs[i] = true
        end
    end
    fault_df = fault_df[fault_idxs, :]
end


function process_faults_from_gis_files(fault_files...;
    layernames=[],
    name="name", dip_dir=:dip_dir, v_ex=:v_ex, e_ex=:e_ex,
    v_rl=:v_rl, e_rl=:e_rl, dip=:dip, hw=:hw, fw=:fw, usd=:usd, lsd=:lsd,
    v_default=0.0, e_default=5.0, usd_default=0.0, lsd_default=15.0, 
    ridge_usd_default=1.0, ridge_lsd_default=2.0, 
    transform_usd_default=1.0, transform_lsd_default=10.0,
    ridge_v_rl=0.0, ridge_e_rl=0.5,
    transform_v_ex=0.0, transform_e_ex=0.5,
    fid=:fid,
    adjust_err_by_dip=false, dip_adj_remainder=0.2,
    fid_drop=[], block_df=:block_df, check_blocks=false,
    subset_in_bounds=false)

    fault_df = load_faults_from_gis_files(fault_files...; layernames)

    dropmissing!(fault_df, :hw)
    dropmissing!(fault_df, :fw)
    fault_df = filter(row -> !(row.hw == ""), fault_df)
    fault_df = filter(row -> !(row.fw == ""), fault_df)

    if subset_in_bounds == true
        fault_df = get_faults_bounding_blocks!(fault_df, block_df)
    end

    faults = process_faults_from_df(fault_df; name=name,
        dip_dir=dip_dir,
        v_ex=v_ex,
        e_ex=e_ex,
        v_rl=v_rl,
        e_rl=e_rl,
        dip=dip,
        hw=hw,
        fw=fw,
        usd=usd,
        lsd=lsd,
        v_default=v_default,
        e_default=e_default,
        usd_default=usd_default,
        lsd_default=lsd_default,
        ridge_usd_default=ridge_usd_default,
        ridge_lsd_default=ridge_lsd_default,
        transform_usd_default=transform_usd_default,
        transform_lsd_default=transform_lsd_default,
        ridge_v_rl=ridge_v_rl,
        ridge_e_rl=ridge_e_rl,
        transform_v_ex=transform_v_ex,
        transform_e_ex=transform_e_ex,
        adjust_err_by_dip=adjust_err_by_dip,
        dip_adj_remainder=dip_adj_remainder,
        fid=fid,
        fid_drop=fid_drop)

    if check_blocks
        faults = Oiler.Utils.check_hw_fw_all(faults, block_df; verbose=true)
        fault_fids = [f.fid for f in faults]
        fault_df = filter(r -> r.fid in fault_fids, fault_df)
    end

    fault_vels = make_vels_from_faults(faults)

    fault_df, faults, fault_vels

end


function get_non_fault_block_bounds(block_df, faults; verbose=false)
    if verbose
        @info "Doing non-fault block boundaries"
        @info "...getting shared block boundary sets"
    end
    shared_bound_sets = Oiler.Utils.get_shared_boundaries(block_df)
    if verbose
        @info "...removing fault segments"
    end
    trimmed_bounds = Oiler.Utils.remove_fault_segs_from_bounds(
        shared_bound_sets, faults
    )

    if verbose
        @info "...sorting segments"
    end
    non_fault_bounds = Oiler.Utils.sort_bound_sets(trimmed_bounds)

    if verbose
        @info "...making boundaries"
    end
    boundaries = []
    for ((fix, mov), traces) in non_fault_bounds
        for (i, trace) in enumerate(traces)
            boundary = Oiler.Boundaries.Boundary(trace, fix, mov, "$fix-$mov-$i")
            push!(boundaries, boundary)
        end
    end
    boundaries
end


function make_vel_from_slip_rate(slip_rate_row, fault_df; err_return_val=1.0,
    weight=1.0, name="name", dip_dir=:dip_dir,
    v_ex=:v_ex, e_ex=:e_ex,
    v_rl=:v_rl, e_rl=:e_rl, dip=:dip, hw=:hw,
    fw=:fw, usd=:usd, lsd=:lsd,
    v_default=0.0, e_default=5.0, usd_default=0.0,
    lsd_default=20.0, fid=:fid)

    fault_fid_type = typeof(fault_df[1, :fid])

    fault_seg = slip_rate_row[:fault_seg]
    # fault_idx = fault_seg# parse(Int, fault_seg)
    if fault_fid_type <: AbstractString
        fault_idx = fault_seg
    else
        fault_idx = parse(fault_fid_type, fault_seg)
    end
    fault_row = @subset(fault_df, :fid .== fault_idx)[1, :]
    fault = row_to_fault(fault_row;
        name=name,
        dip_dir=dip_dir,
        v_ex=v_ex,
        e_ex=e_ex,
        v_rl=v_rl,
        e_rl=e_rl,
        dip=dip,
        hw=hw,
        fw=fw,
        usd=usd,
        lsd=lsd,
        v_default=v_default,
        e_default=e_default,
        usd_default=usd_default,
        lsd_default=lsd_default,
        fid=fid)

    extension_rate = val_nothing_fix(slip_rate_row[:extension_rate])
    extension_err = err_nothing_fix(slip_rate_row[:extension_err];
        return_val=err_return_val, weight=weight)
    dextral_rate = val_nothing_fix(slip_rate_row[:dextral_rate])
    dextral_err = err_nothing_fix(slip_rate_row[:dextral_err];
        return_val=err_return_val, weight=weight)

    ve, vn = Oiler.Faults.fault_slip_rate_to_ve_vn(dextral_rate / weight,
        extension_rate,
        fault.strike)

    ee, en, cen = Oiler.Faults.fault_slip_rate_err_to_ee_en(dextral_err / weight,
        extension_err,
        fault.strike)

    pt = Oiler.IO.get_coords_from_geom(slip_rate_row[:geometry])
    lon = pt[1]
    lat = pt[2]

    VelocityVectorSphere(lon=lon, lat=lat, ve=ve, vn=vn, fix=fault.hw,
        ee=ee, en=en, cen=cen,
        mov=fault.fw, vel_type="geol_slip_rate", name=fault_seg)
end


function make_geol_slip_rate_vels(geol_slip_rate_df, fault_df;
    err_return_val=1.0, weight=1.0, name="name", dip_dir=:dip_dir,
    v_ex=:v_ex, e_ex=:e_ex,
    v_rl=:v_rl, e_rl=:e_rl, dip=:dip, hw=:hw,
    fw=:fw, usd=:usd, lsd=:lsd,
    v_default=0.0, e_default=5.0, usd_default=0.0,
    lsd_default=20.0, fid=:fid,
    warn=false)

    geol_slip_rate_vels = []
    for i in 1:size(geol_slip_rate_df, 1)
        slip_rate_row = geol_slip_rate_df[i, :]
        if (slip_rate_row[:include] == true) | (slip_rate_row[:include] == "1")
            try
                push!(geol_slip_rate_vels, make_vel_from_slip_rate(slip_rate_row,
                    fault_df;
                    err_return_val=err_return_val,
                    weight=weight,
                    name=name,
                    dip_dir=dip_dir,
                    v_ex=v_ex,
                    e_ex=e_ex,
                    v_rl=v_rl,
                    e_rl=e_rl,
                    dip=dip,
                    hw=hw,
                    fw=fw,
                    usd=usd,
                    lsd=lsd,
                    v_default=v_default,
                    e_default=e_default,
                    usd_default=usd_default,
                    lsd_default=lsd_default,
                    fid=fid
                ))
            catch
                warn_msg = "Can't process $slip_rate_row[:fault_seg]"
                if warn
                    @warn warn_msg
                end
            end
        end
    end

    geol_slip_rate_vels = convert(Array{VelocityVectorSphere}, geol_slip_rate_vels)
end


function make_geol_slip_rate_vels!(geol_slip_rate_df, fault_df;
    err_return_val=1.0, weight=1.0, name="name", dip_dir=:dip_dir,
    v_ex=:v_ex, e_ex=:e_ex,
    v_rl=:v_rl, e_rl=:e_rl, dip=:dip, hw=:hw,
    fw=:fw, usd=:usd, lsd=:lsd,
    v_default=0.0, e_default=5.0, usd_default=0.0,
    lsd_default=20.0, fid=:fid,
    warn=false)

    geol_slip_rate_vels = []
    geol_slip_rate_keeps = falses(size(geol_slip_rate_df, 1))
    for i in 1:size(geol_slip_rate_df, 1)
        slip_rate_row = geol_slip_rate_df[i, :]
        if (slip_rate_row[:include] == true) | (slip_rate_row[:include] == "1")
            try
                push!(geol_slip_rate_vels, make_vel_from_slip_rate(slip_rate_row,
                    fault_df;
                    err_return_val=err_return_val,
                    weight=weight,
                    name=name,
                    dip_dir=dip_dir,
                    v_ex=v_ex,
                    e_ex=e_ex,
                    v_rl=v_rl,
                    e_rl=e_rl,
                    dip=dip,
                    hw=hw,
                    fw=fw,
                    usd=usd,
                    lsd=lsd,
                    v_default=v_default,
                    e_default=e_default,
                    usd_default=usd_default,
                    lsd_default=lsd_default,
                    fid=fid
                ))
                geol_slip_rate_keeps[i] = true
            catch
                warn_msg = "Can't process $slip_rate_row[:fault_seg]"
                if warn
                    @warn warn_msg
                end
            end
        end
    end

    geol_slip_rate_df = geol_slip_rate_df[geol_slip_rate_keeps, :]

    for (i, dex_rate) in enumerate(geol_slip_rate_df[!, :dextral_rate])
        if (check_missing(dex_rate))
            if !(check_missing(geol_slip_rate_df[i, :dextral_err]))
                geol_slip_rate_df[i, :dextral_err] = 0.0
            end
        end
    end

    for (i, ex_rate) in enumerate(geol_slip_rate_df[!, :extension_rate])
        if (check_missing(ex_rate))
            if !(check_missing(geol_slip_rate_df[i, :extension_err]))
                geol_slip_rate_df[i, :extension_err] = 0.0
            end
        end
    end

    geol_slip_rate_vels = convert(Array{VelocityVectorSphere}, geol_slip_rate_vels)
    geol_slip_rate_df, geol_slip_rate_vels
end


function get_stoch_slip_rate_dict(stoch_slip_rates)
    n_iters = size(stoch_slip_rates,1)
    rates_out = Dict{String, Dict{String, Array{Float64}}}()
    for flt in stoch_slip_rates[1]
        rates_out[flt.fid] = Dict(
            "dextral_rate"=>Array{Float64}(undef, n_iters),
            "extension_rate"=>Array{Float64}(undef, n_iters),
        )
    end
    for (i, iter) in enumerate(stoch_slip_rates)
        for flt in iter
            rates_out[flt.fid]["dextral_rate"][i] += flt.dextral_rate
            rates_out[flt.fid]["extension_rate"][i] += flt.extension_rate
        end
    end
    rates_out
end


function tri_from_feature(feat)
    props = feat["properties"]
    kps = keys(props)
    
    if "fw" in kps
        fw = string(props["fw"])
    else
        fw = ""
    end

    extra_args = Dict{Symbol, Any}()
    for kw in kps
        if (kw in string.(fieldnames(Oiler.Tris.Tri))) & ((kw != "fw"))
            extra_args[Symbol(kw)] = props[kw]
        end
    end

    tri = Oiler.Tris.Tri(;
        p1=Float64.(feat["geometry"]["coordinates"][1][1]),
        p2=Float64.(feat["geometry"]["coordinates"][1][2]),
        p3=Float64.(feat["geometry"]["coordinates"][1][3]),
        name=string(feat["properties"]["fid"]),
        fw=fw,
        extra_args...
    )
end


function tris_from_geojson(tri_json)
    tris = map(tri_from_feature, tri_json["features"])

    tris = filter(t -> (t.p1 != t.p2) && (t.p1 != t.p3) && (t.p2 != t.p3), tris)

end


function tri_to_feature(tri)

    gj_feature = Dict(
        "type" => "Feature",
        "geometry" => Dict(
            "type" => "Polygon",
            "coordinates" => [[tri.p1, tri.p2, tri.p3, tri.p1]]
        ),
        "properties" => Dict(
            "dip_slip_rate" => round.(tri.dip_slip_rate, digits=3),
            "dip_slip_err" => round.(tri.dip_slip_err, digits=3),
            "strike_slip_rate" => round.(tri.strike_slip_rate, digits=3),
            "strike_slip_err" => round.(tri.strike_slip_err, digits=3),
            "dip_locking_frac" => round.(tri.dip_locking_frac, digits=3),
            "strike_locking_frac" => round.(tri.strike_locking_frac, digits=3),
            "cds" => round.(tri.cds, digits=3),
            "fid" => tri.name,
            "hw" => tri.hw,
            "fw" => tri.fw
        )
    )
end

function tris_to_geojson(tris; name="")
    gj = Dict(
        "type" => "FeatureCollection",
        "name" => name,
        "features" => [tri_to_feature(tri) for tri in tris],
    )
end


function make_tris_from_results(tris, results; default_no_err=0.0)
    tri_rates = results["tri_slip_rates"]
    tri_names = [tri.name for tri in tris]

    function default_return(tri_rate, param; default_val=default_no_err)
        if haskey(tri_rate, param)
            return tri_rate[param]
        else
            return default_val
        end
    end

    tris_out = [Oiler.Tris.Tri(;
        tris[i].p1, tris[i].p2, tris[i].p3,
        dip_slip_rate=tri_rates[name]["dip_slip"],
        strike_slip_rate=tri_rates[name]["strike_slip"],
        dip_slip_err=default_return(tri_rates[name], "dip_slip_err"),
        strike_slip_err=default_return(tri_rates[name], "strike_slip_err"),
        cds=default_return(tri_rates[name], "dsc"),
        dip_locking_frac=tris[i].dip_locking_frac,
        strike_locking_frac=tris[i].strike_locking_frac,
        name=name,
        hw=tris[i].hw,
        fw=tris[i].fw) for (i, name) in enumerate(tri_names)]
end


function write_tri_results_to_gj(tris, results, outfile; name="")
    tris_out = make_tris_from_results(tris, results)
    tri_gj = tris_to_geojson(tris_out; name=name)

    open(outfile, "w") do f
        JSON.print(f, tri_gj)
    end
end


function faults_to_geojson(faults; name="", calc_slip_rate=true, calc_rake=true)
    gj = Dict(
        "type" => "FeatureCollection",
        "name" => name,
        "features" => [fault_to_feature(fault, calc_slip_rate=calc_slip_rate,
            calc_rake=calc_rake) for fault in faults]
    )
end


function fault_to_feature(fault; calc_rake=true, calc_slip_rate=true)
    gj_feature = Dict(
        "type" => "Feature",
        "geometry" => Dict(
            "type" => "LineString",
            "coordinates" => fault.trace'
        ),
        "properties" => Dict(
            "dip" => fault.dip,
            "dip_dir" => fault.dip_dir,
            "extension_rate" => round(fault.extension_rate, digits=3),
            "extension_err" => round(fault.extension_err, digits=3),
            "dextral_rate" => round(fault.dextral_rate, digits=3),
            "dextral_err" => round(fault.dextral_err, digits=3),
            "cde" => round(fault.cde, digits=3),
            "lsd" => fault.lsd,
            "usd" => fault.usd,
            "name" => fault.name,
            "fid" => fault.fid,
            "hw" => fault.hw,
            "fw" => fault.fw
        )
    )

    if calc_slip_rate == true
        slip_rate, slip_rate_err = Oiler.Faults.get_net_slip_rate(fault)
        gj_feature["properties"]["net_slip_rate"] = round(slip_rate, digits=3)
        gj_feature["properties"]["net_slip_rate_err"] = round(slip_rate_err, digits=3)
    end

    if calc_rake == true
        rake, rake_err = Oiler.Faults.get_rake(fault; err=true)
        gj_feature["properties"]["rake"] = round(rake, digits=1)
        gj_feature["properties"]["rake_err"] = round(rake_err, digits=1)
    end

    gj_feature
end


function write_faults_to_gj(faults, outfile; name="", calc_rake=false,
    calc_slip_rate=false)

    fault_gj = faults_to_geojson(faults; name=name, calc_rake=calc_rake,
        calc_slip_rate=calc_slip_rate)

    open(outfile, "w") do f
        JSON.print(f, fault_gj)
    end
end


function write_fault_results_to_gj(results, outfile; name="", calc_rake=false,
    calc_slip_rate=true)
    fault_gj = faults_to_geojson(results["predicted_slip_rates"]; name=name,
        calc_rake=calc_rake, calc_slip_rate=calc_slip_rate)

    open(outfile, "w") do f
        JSON.print(f, fault_gj)
    end
end


function write_stoch_slip_rates_to_gj(stoch_rates, outfile)
    open(outfile, "w") do f
        JSON.print(f, stoch_rates)
    end
end


function write_gnss_vel_results_to_csv(results, vel_groups; name="")
    pred_gnss_df = Oiler.ResultsAnalysis.get_gnss_results(results, vel_groups)
    CSV.write(name, pred_gnss_df)
end


function write_fault_covariance_to_csv(cov_df, outfile)
    CSV.write(outfile, cov_df)
end


function gnss_vel_from_row(row, block; ve=:ve, vn=:vn, ee=:ee, en=:en, cen=:cen,
    fix="1111", name=:station)
    pt = Oiler.IO.get_coords_from_geom(row[:geometry])
    lon = pt[1]
    lat = pt[2]

    vel_field_names = [ve vn ee en cen name]
    vel_fields = [:ve :vn :ee :en :cen :name]
    
    defaults = Dict{Any, Any}(field => 0.0 for field in vel_fields if field != :name)
    defaults[:name] = ""
    vel_vals = Dict()

    for (i, name) in enumerate(vel_field_names)
        vf = vel_fields[i]
        if (name in names(row)) | (string(name) in names(row))
            vel_vals[vf] = row[name]
        else
            vel_vals[vf] = defaults[vf]
        end
    end

    Oiler.VelocityVectorSphere(lon=lon, lat=lat, ve=vel_vals[:ve], vn=vel_vals[:vn], 
        ee=vel_vals[:ee], en=vel_vals[:en], cen=vel_vals[:cen], name=vel_vals[:name],
        fix=fix, mov=string(block), vel_type="GNSS")
end


function make_vels_from_gnss_and_blocks(gnss_df, block_df; ve=:ve, vn=:vn, ee=:ee, en=:en, 
        cen=:cen, fix="1111", name=:station, epsg=102016)

    if epsg != 4326
        block_idx = get_block_idx_for_points(gnss_df, block_df, epsg)
    else
        block_idx = get_block_idx_for_points(gnss_df, block_df)
    end

    gnss_vels = []

    for (i, block) in enumerate(block_idx)
        if !ismissing(block)
            gv = gnss_vel_from_row(gnss_df[i, :], block; ve=ve, vn=vn, ee=ee,
                en=en, cen=cen, fix=fix, name=name)
            push!(gnss_vels, gv)
        end
    end

    gnss_vels = convert(Array{VelocityVectorSphere}, gnss_vels)
end


function load_pred_vels_from_pt_file(pred_vel_pt_file)
    filetype = split(pred_vel_pt_file, ".")[end]
    if filetype in ["geojson", "json"]
        pt_df = gis_vec_file_to_df(pred_vel_pt_file)
    elseif filetype in [".csv"]
        @warn "CSV not yet implemented"
    end

    pt_df
end

function make_vels_for_pred_from_pts(pred_pt_df, block_df; fix="1111", epsg=102016,
    fault_locking=true)

    if fault_locking != true
        @warn "pred vels for no fault locking not implemented yet, here... ignoring"
    end

    n_pts = size(pred_pt_df)[1]

    pred_pt_df[!, "ve"] = zeros(n_pts)
    pred_pt_df[!, "ee"] = zeros(n_pts)
    pred_pt_df[!, "vn"] = zeros(n_pts)
    pred_pt_df[!, "en"] = zeros(n_pts)
    pred_pt_df[!, "station"] = string.(collect(1:n_pts))

    pred_vels = make_vels_from_gnss_and_blocks(pred_pt_df, block_df; epsg=epsg)
end


function make_vels_for_pred_from_pt_file(pred_vel_pt_file, block_df; fix="1111",
    epsg=102016, fault_locking=true)

    pred_pt_df = load_pred_vels_from_pt_file(pred_vel_pt_file)
    pred_vels = make_vels_for_pred_from_pts(pred_pt_df, block_df; fix=fix,
        epsg=epsg, fault_locking=fault_locking)
end


function write_solution_poles(outfile, results, block_df, pole_fix;
    convert_to_sphere=true)
    pole_arr = collect(values(results["poles"]))
    pole_arr = [pole for pole in pole_arr if typeof(pole) == Oiler.PoleCart]

    #rel_poles = [Oiler.Utils.get_path_euler_pole(pole_arr, pole_fix,
    #                                            string(block_df[i, :fid]))
    #            for i in 1:size(block_df, 1)]

    rel_poles = []
    for i in 1:size(block_df, 1)
        try
            pole = Oiler.Utils.get_path_euler_pole(pole_arr, pole_fix,
                string(block_df[i, :fid]))
            push!(rel_poles, pole)
        catch
            mov = string(block_df[i, :fid])
            warn_msg = "can't get pole ($pole_fix, $mov) [$i]"
            @warn warn_msg
        end
    end

    rel_poles = convert(Array{Oiler.PoleCart}, rel_poles)

    CSV.write(outfile, Oiler.IO.poles_to_df(rel_poles,
        convert_to_sphere=convert_to_sphere))
end


function write_results_stats(results, outfile; description=nothing)
    rs = results["stats_info"]

    if !(isnothing(description))
        rs["description"] = description
    end

    open(outfile, "w") do f
        JSON.print(f, rs)
    end
end

function row_to_feature(row; min_dist=0.001, simplify=true,
    check_poly_winding_order=true, epsg=102016)

    if typeof(row[:geometry]) == Oiler.Geom.Polygon
        geom_json = geom_to_geojson(row[:geometry]; simplify=simplify,
            min_dist=min_dist, check_poly_winding_order=check_poly_winding_order,
            epsg=epsg)
    else
        geom_json = geom_to_geojson(row[:geometry])
    end

    feat = Dict{Any,Any}("type" => "Feature",
        "geometry" => geom_json)

    feat["properties"] = Dict()
    for field in names(row)
        if !(field in ["geometry", "geom_type"])
            feat["properties"][field] = row[field]
        end
    end
    feat
end


function geom_to_geojson(geom::Oiler.Geom.Point)
    return Dict("coordinates" => geom.coords,
        "type" => "Point")
end


function geom_to_geojson(geom::Oiler.Geom.LineString)
    return Dict("coordinates":geom.coords,
        "type":"LineString")
end


function geom_to_geojson(geom::Oiler.Geom.Polygon; simplify=false, min_dist=0.0001,
    check_poly_winding_order=true, epsg=102016)
    if simplify
        coords = Oiler.Geom.simplify_polyline(geom.coords, min_dist)
    else
        coords = geom.coords
    end

    coords = Oiler.Utils.matrix_to_rows(coords)

    if check_poly_winding_order
        trans = make_trans_from_wgs84(epsg)
        trans_coords = trans.(coords)

        if Oiler.Geom.check_winding_order(trans_coords) == 1
            reverse!(coords)
        end
    end


    #return Dict("coordinates" => JSON.json([coords]),
    return Dict("coordinates" => [coords],
        "type" => "Polygon")
end


function features_to_geojson(feature_df; name="", min_dist=0.0001, simplify=false,
    check_poly_winding_order=false, epsg=102016)

    features = []


    for i in 1:size(feature_df, 1)
        row = feature_df[i, :]

        # I'm sorry
        if row.fid != "ant"
            push!(features, row_to_feature(row; min_dist=min_dist,
                check_poly_winding_order=check_poly_winding_order,
                epsg=epsg))
        else
            push!(features, row_to_feature(row; min_dist=min_dist,
                check_poly_winding_order=check_poly_winding_order,
                epsg=102019))
        end
    end

    gj = Dict("type" => "FeatureCollection",
        "features" => features)
    if name != ""
        gj["name"] = name
    end

    gj
end


function write_block_df(block_df, outfile; name="", min_dist=0.0001, simplify=false)
    gj = features_to_geojson(block_df; name=name, min_dist=min_dist,
        check_poly_winding_order=true)
    open(outfile, "w") do f
        JSON.print(f, gj)
    end
end


function write_block_centroid_vels_to_csv(results, block_df=nothing; outfile, fix)
    if !(haskey(results, "block_centroids"))
        get_block_centroid_vels(results, block_df; fix=fix)
    end

    CSV.write(outfile, results["block_centroids"])
end


function write_geol_slip_rate_results_to_csv(results; outfile)
    slip_rates = results["data"]["geol_rates_table"]
    CSV.write(outfile, slip_rates)
end


end # module

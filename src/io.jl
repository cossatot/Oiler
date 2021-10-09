module IO

export vel_from_row, load_vels_from_csv, group_vels_by_fix_mov

using ..Oiler
using ..Oiler: VelocityVectorSphere, PoleSphere, PoleCart, pole_cart_to_sphere


using Logging
import Base.Threads.@threads

using CSV
using JSON
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


function gis_vec_file_to_df(filename::AbstractString; 
                            layername::AbstractString="")
    dataset = AG.read(filename)
    if layername == ""
        layer = AG.getlayer(dataset, 0)
    else
        layer = AG.getlayer(dataset, layername)
    end

    dataframe = DataFrame(layer)
    
    if !("geometry" in names(dataframe))
        rename!(dataframe, "" => :geometry)
    end

    if !("fid" in names(dataframe))
        nfeat = AG.nfeature(layer)
        fids = []
        for nf in 1:nfeat
            feature = AG.unsafe_nextfeature(layer)
            push!(fids, AG.getfid(feature))
        end
        dataframe[!, "fid"] = fids
    end

    dataframe
end


function gis_vec_file_to_df_old(filename::AbstractString; layername="")
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


# function get_block_idx_for_point_list(point_list, block_df; epsg=4326)
#   point_geoms = map(p->AG.createpoint(p[1],p[2]), point_list)
#
#   if epsg != 4326
#       point_geoms = map(x->AG.reproject(x,
#                           ProjString("+proj=longlat +datum=WGS84 +no_defs"), 
#                           EPSG(epsg)),
#                   point_geoms)
#   end
#
#   idxs = Array{Any}(missing, length(point_geoms))
#
#   for i_b in 1:size(block_df, 1)
#       block_geom = AG.clone(block_df[i_b, :geometry])
#       if epsg != 4326
#           block_geom = AG.reproject(block_geom, 
#               ProjString("+proj=longlat +datum=WGS84 +no_defs"), EPSG(epsg))
#       end
#
#       block_fid = block_df[i_b, :fid]
#
#       for (i_p, pg) in point_geoms
#           if AG.contains(block_geom, pg)
#               idxs[i_p] = block_fid
#           end
#       end
#   end
#   idxs
# end


function get_block_idx_for_point(point, block_df; epsg=4326)
    point_geom = AG.createpoint(point[1], point[2])

    if epsg != 4326
        point_geom = AG.reproject(point_geom, 
            ProjString("+proj=longlat +datum=WGS84 +no_defs"), EPSG(epsg))
    end

    idx = missing

    for i_b in 1:size(block_df, 1)
        block_geom = AG.clone(block_df[i_b, :geometry])
        if epsg != 4326
            block_geom = AG.reproject(block_geom, 
                ProjString("+proj=longlat +datum=WGS84 +no_defs"), EPSG(epsg))
        end

        block_fid = block_df[i_b, :fid]

        if AG.contains(block_geom, point_geom)
            idx = block_fid
        end
    end
    idx
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


function row_to_fault(row; name="name", dip_dir=:dip_dir, v_ex=:v_ex, e_ex=:e_ex,
        v_rl=:v_rl, e_rl=:e_rl, dip=:dip, hw=:hw, fw=:fw, usd=:usd, lsd=:lsd,
        v_default=0., e_default=5., usd_default=0., lsd_default=20., fid=:fid)

    trace = Oiler.IO.get_coords_from_geom(row[:geometry])
    if name in names(row)
        _name = row[:name]
    else
        _name = ""
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
        fid=row[fid]
        )
end


function process_faults_from_df(fault_df; name="name", dip_dir=:dip_dir, v_ex=:v_ex, 
                                e_ex=:e_ex, v_rl=:v_rl, e_rl=:e_rl, dip=:dip,
                                hw=:hw, fw=:fw, usd=:usd, lsd=:lsd,
                                v_default=0., e_default=5., usd_default=0.,
                                lsd_default=20., fid=:fid, fid_drop=:fid_drop)

    faults = []
    for i in 1:size(fault_df, 1)
        # try
        push!(faults, Oiler.IO.row_to_fault(fault_df[i,:]; 
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
                                                fid=fid) ) 
        # catch
        #    warn_msg = "Problem with $fault_df[i,:]"
        #    @warn warn_msg 
        # end 
    end
    faults = [f for f in faults if f.fw != ""]
    faults = [f for f in faults if f.hw != ""]
    faults = [f for f in faults if f.fid != fid_drop]

    faults
end


function make_vels_from_faults(faults)
    fault_vels = reduce(vcat, map(Oiler.fault_to_vels, faults))
    fault_vels
end



function load_faults_from_gis_files(fault_files... ; layernames=[])
    
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


function get_blocks_in_bounds!(block_df, bound_df; epsg=0)
    bound_orig = bound_df[1,:geometry]
    block_geoms_orig = block_df[:,:geometry]

    if epsg != 0
        bound = AG.reproject(AG.clone(bound_orig), 
            ProjString("+proj=longlat +datum=WGS84 +no_defs"), EPSG(epsg))
        block_geoms = [AG.reproject(AG.clone(g), 
                                    ProjString("+proj=longlat +datum=WGS84 +no_defs"), 
                                    EPSG(epsg))
                       for g in block_geoms_orig]
    else
        bound = bound_orig
        block_geoms = block_geoms_orig
    end

    block_idxs = falses(length(block_geoms))
    for (i, block_geom) in enumerate(block_geoms)
        if AG.intersects(bound, block_geom)
            block_idxs[i] = true
        end
    end
    block_df[block_idxs,:]
end


function get_blocks_and_faults_in_bounds!(block_df, fault_df, bound_df)
    block_df = get_blocks_in_bounds!(block_df, bound_df)
    fault_df = get_faults_bounding_blocks!(fault_df, block_df)

    block_df, fault_df
end


function get_faults_bounding_blocks!(fault_df, block_df)
    block_fids = string.(block_df[:,:fid])
    fault_idxs = falses(size(fault_df, 1))
    for i in 1:length(fault_idxs)
        if (fault_df[i, :hw] in block_fids) & (fault_df[i,:fw] in block_fids)
            fault_idxs[i] = true
        end
    end
    fault_df = fault_df[fault_idxs,:]
end


function process_faults_from_gis_files(fault_files... ; 
        layernames=[],
        name="name", dip_dir=:dip_dir, v_ex=:v_ex, e_ex=:e_ex,
        v_rl=:v_rl, e_rl=:e_rl, dip=:dip, hw=:hw, fw=:fw, usd=:usd, lsd=:lsd,
        v_default=0., e_default=5., usd_default=0., lsd_default=20., fid=:fid,
        fid_drop=:fid_drop, block_df=:block_df,
        subset_in_bounds=false)
    
    fault_df = load_faults_from_gis_files(fault_files ...; layernames)
    
    if subset_in_bounds == true
        fault_df = get_faults_bounding_blocks!(fault_df, block_df)
    end

    dropmissing!(fault_df, :hw)
    dropmissing!(fault_df, :fw)

    fault_df = filter(row -> !(row.hw == ""), fault_df)
    fault_df = filter(row -> !(row.fw == ""), fault_df)

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
                                    fid=fid,
                                    fid_drop=fid_drop)

    fault_vels = make_vels_from_faults(faults)

    fault_df, faults, fault_vels

end


function make_vel_from_slip_rate(slip_rate_row, fault_df; err_return_val=1., 
                                 weight=1., name="name", dip_dir=:dip_dir, 
                                 v_ex=:v_ex, e_ex=:e_ex,
                                 v_rl=:v_rl, e_rl=:e_rl, dip=:dip, hw=:hw, 
                                 fw=:fw, usd=:usd, lsd=:lsd,
                                 v_default=0., e_default=5., usd_default=0., 
                                 lsd_default=20., fid=:fid)
    
    fault_fid_type = typeof(fault_df[1,:fid])
    
    fault_seg = slip_rate_row[:fault_seg]
    # fault_idx = fault_seg# parse(Int, fault_seg)
    if fault_fid_type <: AbstractString
        fault_idx = fault_seg
    else
        fault_idx = parse(fault_fid_type, fault_seg)
    end
    fault_row = @where(fault_df, :fid .== fault_idx)[1,:]
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
                         mov=fault.fw, vel_type="fault", name=fault_seg) 
end


function make_geol_slip_rate_vels(geol_slip_rate_df, fault_df;
        err_return_val=1., weight=1., name="name", dip_dir=:dip_dir, 
        v_ex=:v_ex, e_ex=:e_ex,
        v_rl=:v_rl, e_rl=:e_rl, dip=:dip, hw=:hw, 
        fw=:fw, usd=:usd, lsd=:lsd,
        v_default=0., e_default=5., usd_default=0., 
        lsd_default=20., fid=:fid,
        warn=false)

    geol_slip_rate_vels = []
    for i in 1:size(geol_slip_rate_df, 1)
        slip_rate_row = geol_slip_rate_df[i,:]
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
        err_return_val=1., weight=1., name="name", dip_dir=:dip_dir, 
        v_ex=:v_ex, e_ex=:e_ex,
        v_rl=:v_rl, e_rl=:e_rl, dip=:dip, hw=:hw, 
        fw=:fw, usd=:usd, lsd=:lsd,
        v_default=0., e_default=5., usd_default=0., 
        lsd_default=20., fid=:fid,
        warn=false)

    geol_slip_rate_vels = []
    geol_slip_rate_keeps = falses(size(geol_slip_rate_df, 1))
    for i in 1:size(geol_slip_rate_df, 1)
        slip_rate_row = geol_slip_rate_df[i,:]
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

    geol_slip_rate_df = geol_slip_rate_df[geol_slip_rate_keeps,:]
    geol_slip_rate_vels = convert(Array{VelocityVectorSphere}, geol_slip_rate_vels)
    geol_slip_rate_df, geol_slip_rate_vels
end




function tri_from_feature(feat)
    props = feat["properties"]
    kps = keys(props)
    if "fw" in kps
        fw = string(props["fw"])
    else fw = ""
    end

    tri = Oiler.Tris.Tri(;
       p1=Float64.(feat["geometry"]["coordinates"][1][1]),
       p2=Float64.(feat["geometry"]["coordinates"][1][2]),
       p3=Float64.(feat["geometry"]["coordinates"][1][3]),
       name=string(feat["properties"]["fid"]),
       fw=fw
       )
end


function tris_from_geojson(tri_json)
    tris = map(tri_from_feature, tri_json["features"])
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


function make_tris_from_results(tris, results)
    tri_rates = results["tri_slip_rates"]
    tri_names = [tri.name for tri in tris]
    tris_out = [Oiler.Tris.Tri(;
                tris[i].p1, tris[i].p2, tris[i].p3,
                dip_slip_rate=tri_rates[name]["dip_slip"],
                strike_slip_rate=tri_rates[name]["strike_slip"],
                dip_slip_err=tri_rates[name]["dip_slip_err"],
                strike_slip_err=tri_rates[name]["strike_slip_err"],
                cds=tri_rates[name]["dsc"],
                dip_locking_frac=tris[i].dip_locking_frac,
                strike_locking_frac=tris[i].strike_locking_frac,
                name=name,
                hw=tris[i].hw,
                fw=tris[i].fw,

                ) for (i, name) in enumerate(tri_names)]
end


function write_tri_results_to_gj(tris, results, outfile; name="")
    tris_out = make_tris_from_results(tris, results)
    tri_gj = tris_to_geojson(tris_out; name=name)

    open(outfile, "w") do f
        JSON.print(f, tri_gj)
    end
end


function faults_to_geojson(faults; name="")
    gj = Dict(
        "type" => "FeatureCollection",
        "name" => name,
        "features" => [fault_to_feature(fault) for fault in faults]
    )
end


function fault_to_feature(fault)
    gj_feature = Dict(
        "type" => "Feature",
        "geometry" => Dict(
            "type" => "LineString",
            "coordinates" => fault.trace'
        ),
        "properties" => Dict(
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
end


function write_fault_results_to_gj(results, outfile; name="")
    fault_gj = faults_to_geojson(results["predicted_slip_rates"]; name=name)
    
    open(outfile, "w") do f
        JSON.print(f, fault_gj)
    end
end


function write_gnss_vel_results_to_csv(results, vel_groups; name="")
    pred_gnss_df = Oiler.ResultsAnalysis.get_gnss_results(results, vel_groups)
    CSV.write(name, pred_gnss_df)
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


function write_solution_poles(outfile, results, block_df, pole_fix;
                              convert_to_sphere=true)
    pole_arr = collect(values(results["poles"]))
    pole_arr = [pole for pole in pole_arr if typeof(pole) == Oiler.PoleCart]

    rel_poles = [Oiler.Utils.get_path_euler_pole(pole_arr, pole_fix,
                                                string(block_df[i, :fid]))
                for i in 1:size(block_df, 1)]

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


function row_to_feature(row; coord_digits=5)
    feat = JSON.parse(AG.toJSON(row[:geometry], 
                                COORDINATE_PRECISION=coord_digits))
    feat["properties"] = Dict()
    for field in names(feat)
        if field != "geometry"
            feat["properties"][field] = row[field]
        end
    end
    feat
end


function features_to_geojson(feature_df; name="", coord_digits=5)
    gj = Dict("type" => "FeatureCollection", 
              features => [row_to_feature(feature_df[i,:]; coord_digits=coord_digits)
                           for i in 1:size(feature_df, 1)])
    if name != ""
        gj["name"] = name
    end
    gj
end


function write_block_df(block_df, outfile; name="", coord_digits=5)
    gj = features_to_geojson(block_df; name=name, coord_digits=coord_digits)
    open(outfile, "w") do f
        JSON.print(f, gj)
    end


end # module
module Utils

export predict_vels_from_poles, find_vel_cycles, diagonalize_matrices,
    get_gnss_vels, get_coords_from_vel_array

using ..Oiler
using ..Oiler: VelocityVectorSphere, PoleCart, PoleSphere, build_PvGb_from_vels,
    build_vel_column_from_vels, add_poles, pole_sphere_to_cart

using Logging
using SparseArrays
using LinearAlgebra
import Base.Threads.@threads

using ThreadsX
using Setfield
using DataFrames


function diagonalize_matrices(matrices; sparse_output=true)

    rowz = [size(m, 1) for m in matrices]
    colz = [size(m, 2) for m in matrices]

    n_rows = sum(rowz)
    n_cols = sum(colz)

    if sparse_output == true
        big_mat = spzeros(n_rows, n_cols)
    else
        big_mat = zeros(n_rows, n_cols)
    end

    i_row = 1
    i_col = 1

    for (i, mat) in enumerate(matrices)
        row_end = i_row + rowz[i] - 1
        col_end = i_col + colz[i] - 1

        big_mat[i_row:row_end, i_col:col_end] = mat

        i_row = row_end + 1
        i_col = col_end + 1
    end
    return big_mat
end


function sparse_to_dict(sparray::SparseMatrixCSC)
    # (rows, cols, vals) = findnz(sparray)
    nz_inds = Tuple(findall(!iszero, sparray))
    rows = [nzi[1] for nzi in nz_inds]
    cols = [nzi[2] for nzi in nz_inds]
    Dict((r, cols[i]) => sparray[r, cols[i]] for (i, r) in enumerate(rows))
end


function dict_to_sparse(sparse_dict::Dict)
    rowz = [] # Array{Int}(length(sparse_dict))
    colz = [] # Array{Int}(length(sparse_dict))
    valz = [] # Array{}(length(sparse_dict))
    for ((r, c), v) in sparse_dict
        push!(rowz, r)
        push!(colz, c)
        push!(valz, v)
    end
    sparse(rowz, colz, valz)
end


function matrix_to_rows(matrix::Matrix)
    return [matrix[i, :] for i in 1:size(matrix, 1)]
end



"""
    lin_indep_cols(X, tol)

    Finds and returns linearly independent columns of X.

    # Arguments:
    - X: Matrix with potentially linearly dependent columns.
    - tol: Rank estimation tolerance.

    Based on code by Matt J: 
    https://www.mathworks.com/matlabcentral/answers/
    108835-how-to-get-only-linearly-independent-rows-
    in-a-matrix-or-to-remove-linear-dependency-b-w-rows-in-a-m
"""
function lin_indep_cols(X; tol=1e-10)
    if ~(true) # supposed to check for all zeros here
        idx = []
    else
        Xm = Matrix(X)

        F = qr(Xm, Val(true))
        Q = F.Q
        R = F.R
        P = F.p

        if !(1 in size(R))
            diagr = abs.(diag(R))
        else
            diagr = R[1]
        end

        r = findlast(diagr .>= tol * diagr[1])
        idx = sort(P[1:r])
    end
    idx
end


function lin_indep_rows(X; tol=1e-10)
    lin_indep_cols(X'; tol=tol)
end


function replacenan!(A::Array{Float64}; replacement::Float64=0.)
    @info "\tdense"
    for i in eachindex(A)[1]
        @inbounds A[i] = ifelse(isnan(A[i]), replacement, A[i])
    end
end

function replacenan!(A::SparseMatrixCSC{Float64}; replacement::Float64=0.)
    @info "\tsparse"
    for i in findnz(A)[1]
        @inbounds A[i] = ifelse(isnan(A[i]), replacement, A[i])
    end
end

function replacenan_multi!(A::Array{Float64}; replacement::Float64=0.)
    @threads for i in eachindex(A)
        @inbounds A[i] = ifelse(isnan(A[i]), replacement, A[i])
    end
end

function make_digraph_from_vels(vels::VelocityVectorSphere...)
    vel_array = collect(vels)
    make_digraph_from_vels(vel_array)
end


"""
    make_digraph_from_vels(vels)

Makes a directed graph from a group of velocity vectors. Each velocity vector
(of type VelocityVectorSphere) has 'fix' and 'mov' fields that describe the
velocity of a point or block (mov) relative to another point or block (fix). The
directed graph maps out these relationships to determine constraints on block
motions, in particular ensuring closure of velocity circuits. Edges are not
duplicated, so if multiple velocities constrain relative motion between two
blocks, only one edge is given in the resulting graph.

The graph is in an adjacency list format.

# Arguments

# Returns


"""
function make_digraph_from_vels(vels::Array{VelocityVectorSphere})
    vel_graph = Dict{String,Array{String}}()

    for vel in vels
        if vel.fix == vel.mov
            # It's unclear to me how to handle velocities that have the
            # same fix and mov.  This is OK for GNSS velocities but not for
            # faults, in principle. But I do not know how it will affect the
            # results.
            if vel.vel_type == "GNSS"
                if !(vel.mov in vel_graph[vel.fix])
                    push!(vel_graph[vel.fix], vel.mov)
                end
            else
                println(vel.vel_type)
                @warn "$vel has same fix and mov, leaving it out."
            end
        elseif haskey(vel_graph, vel.fix)
            if !(vel.mov in vel_graph[vel.fix])
                push!(vel_graph[vel.fix], vel.mov)
            end
        else
            vel_graph[vel.fix] = [vel.mov]
        end
    end
    vel_graph
end


function make_digraph_from_poles(poles::Union{Array{PoleCart},Array{PoleSphere}})
    pole_graph = Dict{String,Array{String}}()

    for pole in poles
        if haskey(pole_graph, pole.fix)
            if !(pole.mov in pole_graph[pole.fix])
                push!(pole_graph[pole.fix], pole.mov)
            end
        else
            pole_graph[pole.fix] = [pole.mov]
        end
    end
    pole_graph
end


function make_digraph_from_tuples(tups)
    graph = Dict()

    for (t1, t2) in tups
        if haskey(graph, t1)
            if !(t2 in graph[t1])
                push!(graph[t1], t2)
            end
        else
            graph[t1] = [t2]
        end
    end
    graph
end


function make_ugraph_from_digraph(digraph::Dict)
    ug = Dict{String,Array{String}}()

    for (fix, movs) in digraph
        for mov in movs
            if !haskey(ug, fix)
                ug[fix] = [mov]
            elseif !(mov in ug[fix])
                push!(ug[fix], mov)
            end
            if !haskey(ug, mov)
                ug[mov] = [fix]
            elseif !(fix in ug[mov])
                push!(ug[mov], fix)
            end
        end
    end
    ug
end


function find_tricycles(graph::Dict; reduce::Bool=true)

    all_tris = []

    for (from, tos) in graph
        for to in tos
            try
                if to != from
                    for next in graph[to]
                        if next != from
                            for fin in graph[next]
                                if fin == from
                                    path = [from, to, next]
                                    push!(all_tris, path)
                                end
                            end
                        end
                    end
                end
            catch e
                if isa(e, KeyError)
                    err_msg = "problem with ($from, $to)"
                    @warn err_msg
                end
            end
        end
    end

    if reduce == true
        sets = []
        unique_tris = []

        for tri in all_tris
            if !(Set(tri) in sets)
                push!(sets, Set(tri))
                push!(unique_tris, tri)
            end
        end
        return unique_tris
    else
        return all_tris
    end
end


function get_cycle_inds(vel_group_keys, cycle_tup)

    v_ind = findall(x -> x == cycle_tup, vel_group_keys)
    if v_ind != []
        res = Dict("ind" => v_ind[1], "val" => 1.0)
    else
        cycle_tup = Base.reverse(cycle_tup)
        v_ind = findall(x -> x == (cycle_tup), vel_group_keys)
        if v_ind != []
            res = Dict("ind" => v_ind[1], "val" => -1.0)
        else
            println("v apparently not in cycle")
        end
    end
    res
end


function find_vel_cycles(vel_group_keys)

    cycle_inds = Dict()
    dg = make_digraph_from_tuples(vel_group_keys)
    tricycles = find_tricycles(make_ugraph_from_digraph(dg))

    for (i, tri) in enumerate(tricycles)
        cycle_inds[i] = Dict()
        for ctup in [(tri[1], tri[2]), (tri[2], tri[3]), (tri[3], tri[1])]
            cycle_inds[i][ctup] = get_cycle_inds(vel_group_keys, ctup)
        end
    end
    cycle_inds
end


"""
    flat(x)
Flattens nested arrays of arrays
"""
flat(x, y=vcat(x...)) = x == y ? x : flat(y)


"""
    find_shortest_path(graph, start, stop)

Finds the shortest path between vertices `start` and `stop` given the
edges in the (directed, unweighted) `graph` using a breadth-first-search
strategy.

# Arguments
- `graph`: A dictionary with keys for each vertex, and values of arrays of
  all the vertices linked to the key vertex.

- `start`: Vertex to start from.

- `stop`: Vertex to stop on.

# Returns
Array of vertices representing the shortest path in order from `start` to 
`stop`.
"""
function find_shortest_path(graph::Dict{String,Array{String}},
    start::String, stop::String)

    dist = Dict{String,Any}(start => [start])
    q = [start]
    while length(q) > 0
        at = popfirst!(q)
        for next in graph[at]
            if !haskey(dist, next)
                dist[next] = [dist[at], next]
                push!(q, next)
            end
        end
    end
    flat(dist[stop])
end


function get_pole_path(poles::Array{PoleCart}, path::Array{String})

    steps = length(path) - 1
    pole_path = Array{PoleCart,1}(undef, steps)

    for i = 1:steps
        place = path[i]
        next = path[i+1]

        for pole in poles
            if pole.fix == place && pole.mov == next
                pole_path[i] = pole
            elseif pole.fix == next && pole.mov == place
                pole_path[i] = -pole
            end
        end
    end
    pole_path
end


function get_path_euler_pole(poles::Array{PoleCart,1}, fix::String,
    mov::String)

    if fix == mov
        final_pole = PoleCart(x=0.0, y=0.0, z=0.0, fix=fix, mov=mov)

    elseif length([p for p in poles if (p.fix == fix) & (p.mov == mov)]) == 1
        final_pole = [p for p in poles if (p.fix == fix) & (p.mov == mov)][1]

    elseif length([p for p in poles if (p.fix == mov) & (p.mov == fix)]) == 1
        final_pole = [p for p in poles if (p.fix == fix) & (p.mov == mov)][1]
        final_pole = -final_pole

    else
        vel_dg = make_digraph_from_poles(poles)
        vel_ug = make_ugraph_from_digraph(vel_dg)

        path = find_shortest_path(vel_ug, fix, mov)

        pole_path = get_pole_path(poles, path)

        final_pole = add_poles(pole_path)
    end
    final_pole
end


function get_path_euler_pole(poles::Dict, fix::String, mov::String)

    if haskey(poles, (fix, mov))
        final_pole = poles[(fix, mov)]
    elseif haskey(poles, (mov, fix))
        final_pole = -poles[(mov, fix)]
    else
        pole_arr = [v for v in values(poles)]
        final_pole = get_path_euler_pole(pole_arr, fix, mov)
    end
    final_pole
end


function predict_vels_from_poles(block_things::Dict{String,AbstractArray},
    poles::Array{PoleCart,1})

    @warn "predict_vels_from_poles  has no error propagation"
    pole_list = [get_path_euler_pole(poles, fix, mov) for (fix, mov) in
                 block_things["keys"]]

    pole_vec = reduce(vcat, [[p.x; p.y; p.z] for p in pole_list])

    V_pred = block_things["PvGb"] * pole_vec

    Vn_pred = V_pred[1:3:end]
    Ve_pred = V_pred[2:3:end]
    # Vu_pred = V_pred[3:3:end]

    (Ve_pred, Vn_pred)
end

function predict_vels_from_poles(block_things::Dict{String,AbstractArray},
    poles::Array{PoleSphere})

    cpoles = [pole_sphere_to_cart(pole) for pole in poles]

    predict_vels_from_poles(block_things, cpoles)
end


function get_vel_vec_from_pole(vel::VelocityVectorSphere, pole::PoleCart)
    PvGb = Oiler.BlockRotations.build_PvGb_vel(vel)
    vel_vec = PvGb * [pole.x; pole.y; pole.z]
end


function get_vel_vec_from_pole(vel::VelocityVectorSphere, pole::PoleSphere)
    get_vel_vec_from_pole(vel, Oiler.BlockRotations.pole_sphere_to_cart(pole))
end


function check_vel_closures(poles; tol=1e-5)
    v_keys = sort(collect(Tuple(keys(poles))))
    block_digraph = make_digraph_from_tuples(v_keys)
    block_ugraph = make_ugraph_from_digraph(block_digraph)
    cycles = find_tricycles(block_ugraph)

    all_good = true

    for cycle in cycles
        pole_12 = get_path_euler_pole(poles, cycle[1], cycle[2])
        pole_23 = get_path_euler_pole(poles, cycle[2], cycle[3])
        pole_31 = get_path_euler_pole(poles, cycle[3], cycle[1])

        p_sum = Oiler.pole_cart_to_sphere(pole_12 + pole_23 + pole_31)

        if p_sum.rotrate > tol
            rr = p_sum.rotrate
            warn_msg = "$cycle does not close: $rr"
            all_good = false
            @warn warn_msg
        end
    end

    if all_good == true
        @info "All poles close"
    end
    all_good
end


function get_fault_slip_rates_from_poles(faults, poles; use_path=true)
    rates = []
    for fault in faults
        fault_key = (fault.fw, fault.hw)
        if haskey(poles, fault_key)
            pole = poles[fault_key]
            push!(rates, Oiler.Faults.get_fault_slip_rate_from_pole(fault, pole))
        elseif haskey(poles, reverse(fault_key))
            pole = poles[reverse(fault_key)]
            push!(rates, Oiler.Faults.get_fault_slip_rate_from_pole(fault, pole))
        else
            if use_path
                try
                    pole = get_path_euler_pole(poles, fault.fw, fault.hw)
                    push!(rates, Oiler.Faults.get_fault_slip_rate_from_pole(
                        fault, pole))
                catch
                    @warn "Can't make pole for $fault_key"
                    push!(rates, (NaN, NaN))
                end
            else
                @warn "No pole found for $fault_key"
                push!(rates, (NaN, NaN))
            end
        end
    end
    rates
end


function group_faults(faults, vel_group_keys)
    fault_groups = Dict(k => [] for k in vel_group_keys)

    for fault in faults
        for k in vel_group_keys
            if fault.fw in k && fault.hw in k
                push!(fault_groups[k], fault)
            end
        end
    end

    fault_groups = Dict(k => v for (k, v) in fault_groups if v != [])
end


function get_n_blocks_from_vel_groups(vel_groups)
    vgk = collect(keys(vel_groups))

    vel_refs = []
    for vg in vgk
        r1 = vg[1]
        r2 = vg[2]

        if !(r1 in vel_refs)
            push!(vel_refs, r1)
        end
        if !(r2 in vel_refs)
            push!(vel_refs, r2)
        end
    end

    n_blocks = length(vel_refs) - 1 # assuming 1 reference frame
end



"""
    get_fault_vels(vel_groups)

Collects all of the faults from the velocity groups, with
some ancillary metadata.

# Arguments

# Returns

"""
function get_fault_vels(vel_groups)
    faults = []
    row_set_num = 0
    vg_keys = sort(collect(Tuple(keys(vel_groups))))

    for (i, key) in enumerate(vg_keys)
        group = vel_groups[key]
        col_idx = 3 * (i - 1) + 1
        for vel in group
            row_set_num += 1
            if vel.vel_type == "fault"
                row_idx = 3 * (row_set_num - 1) + 1
                vel_idx = [row_idx:row_idx+2, col_idx:col_idx+2]
                vd = Dict()
                vd["fault"] = vel
                vd["idx"] = vel_idx
                push!(faults, vd)
            end # if
        end # for
    end # for
    faults
end


"""
    get_gnss_vels(vel_groups)

Collects all of the GNSS-derived velocities from the velocity groups, with
some ancillary metadata.

# Arguments

# Returns

"""
function get_gnss_vels(vel_groups)
    gnss_vels = []
    row_set_num = 0
    vg_keys = sort(collect(Tuple(keys(vel_groups))))

    for (i, key) in enumerate(vg_keys)
        group = vel_groups[key]
        col_idx = 3 * (i - 1) + 1
        for vel in group
            row_set_num += 1
            if vel.vel_type == "GNSS"
                row_idx = 3 * (row_set_num - 1) + 1
                vel_idx = [row_idx:row_idx+2, col_idx:col_idx+2]
                vd = Dict()
                vd["vel"] = vel
                vd["idx"] = vel_idx
                push!(gnss_vels, vd)
            end # if
        end # for
    end # for
    gnss_vels
end


function get_geol_slip_rate_vels(vel_groups)
    geol_slip_rate_vels = []
    row_set_num = 0
    vg_keys = sort(collect(Tuple(keys(vel_groups))))

    for (i, key) in enumerate(vg_keys)
        group = vel_groups[key]
        col_idx = 3 * (i - 1) + 1
        for vel in group
            row_set_num += 1
            if vel.vel_type == "geol_slip_rate"
                row_idx = 3 * (row_set_num - 1) + 1
                vel_idx = [row_idx:row_idx+2, col_idx:col_idx+2]
                vd = Dict()
                vd["vel"] = vel
                vd["idx"] = vel_idx
                push!(geol_slip_rate_vels, vd)
            end # if
        end # for
    end # for
    geol_slip_rate_vels
end


"""
    get_coords_from_vel_array(vels)

Returns (lons, lats) as arrays of floats from an array of `VelocityVectorSphere`
"""
function get_coords_from_vel_array(vels)
    lats = [v.lat for v in vels]
    lons = [v.lon for v in vels]

    (lons, lats)
end


function make_df_from_vel_array(vels)
    lons, lats = get_coords_from_vel_array(vels)
    vel_df = DataFrame(lon=lons, lat=lats)

    vel_df.name = [v.name for v in vels]
    vel_df.fix = [v.fix for v in vels]
    vel_df.mov = [v.mov for v in vels]
    vel_df.ve = [v.ve for v in vels]
    vel_df.vn = [v.vn for v in vels]
    vel_df.ee = [v.ee for v in vels]
    vel_df.en = [v.en for v in vels]

    vel_df
end


function make_gnss_df_from_vel_groups(vel_groups)
    vel__ = get_gnss_vels(vel_groups)
    vels = [v["vel"] for v in vel__]
    make_df_from_vel_array(vels)
end


function tri_priors_from_pole(tris, pole; locking_fraction=1.0,
    err_coeff=1.0, depth_adjust=false, depth_max=80.0)

    err_coeff = 1.0 / err_coeff

    #tri_rates = ThreadsX.map(x -> Oiler.Tris.get_tri_rate_from_pole(x, pole), tris)
    #tri_rates = map(x -> Oiler.Tris.get_tri_rate_from_pole(x, pole), tris)
    tri_rates = []
    for (i, tri) in enumerate(tris)
        push!(tri_rates, Oiler.Tris.get_tri_rate_from_pole(tri, pole))
    end

    function set_tri_rates(tri, rate_array, locking_fraction, err_coeff, depth_adjust, depth_max)

        if depth_adjust
            depth = -Oiler.Tris.get_tri_center(tri)[3]
            if depth < depth_max
                depth_fraction = (depth_max - depth) / depth_max
            else
                depth_fraction = 0.0
            end
        else
            depth_fraction = 1.0
        end

        tri = @set tri.dip_slip_rate = rate_array[1] * locking_fraction * depth_fraction
        tri = @set tri.dip_slip_err = max(rate_array[2], 1.0) * err_coeff
        tri = @set tri.strike_slip_rate = rate_array[3] * locking_fraction * depth_fraction
        tri = @set tri.strike_slip_err = max(rate_array[4], 1.0) * err_coeff
        tri = @set tri.cds = rate_array[5]
    end

    @threads for (i, tri) in collect(enumerate(tris))
        tris[i] = set_tri_rates(tri, tri_rates[i], locking_fraction, err_coeff,
            depth_adjust, depth_max)
    end

    tris
end


function find_first_nz_per_row(matrix)
    first_idxs = []
    nr, nc = size(matrix)
    for i in 1:nr
        for j in 1:nc
            if matrix[i, j] != 0.0
                push!(first_idxs, (i, j))
                break
            end
        end
    end
    first_idxs
end


function sort_sparse_matrix(matrix)
    first_idxs = find_first_nz_per_row(matrix)

    first_idxs = sort!(first_idxs, by=x -> x[2])
    rows = [x[1] for x in first_idxs]

    matrix[rows, :]
end


function get_block_centroids(block_geoms; epsg=102016)
    centroids = map(x -> Oiler.Geom.get_polygon_centroid(x; epsg=epsg), block_geoms)
end


function get_block_centroids(block_df::DataFrame; epsg=102016)
    centroids = get_block_centroids(block_df.geometry; epsg=epsg)
end


function get_block_centroid_velocities(block_df, poles; fix, epsg=102016)
    # returns velocties at block centroid without locking effects
    centroids = map(x -> Oiler.Geom.get_polygon_centroid(x; epsg=epsg), block_df.geometry)

    rows_keep = []
    vels = []

    for (i, centroid) in enumerate(centroids)
        #try
        lon, lat = centroid.coords[1], centroid.coords[2]
        pole = Oiler.Utils.get_path_euler_pole(poles, fix, string(block_df[i, :fid]))
        vel = Oiler.BlockRotations.predict_block_vel(lon, lat, pole)
        push!(vels, vel)
        push!(rows_keep, i)
        #catch
        #end
    end

    fids = string.([block_df[i, :fid] for i in rows_keep])
    df = DataFrame("fid" => fids)
    #df[!, "geometry"] = [centroids[i] for i in rows_keep]
    df[!, "lon"] = [centroids[i].coords[1] for i in rows_keep]
    df[!, "lat"] = [centroids[i].coords[2] for i in rows_keep]
    df[!, "ve"] = [vel.ve for vel in vels]
    df[!, "vn"] = [vel.vn for vel in vels]
    df[!, "ee"] = [vel.ee for vel in vels]
    df[!, "en"] = [vel.en for vel in vels]
    df[!, "cen"] = [vel.cen for vel in vels]

    df
end


function line_to_segs(coords; digits=4)
    line = matrix_to_rows(coords)
    [round.([line[i][1] line[i][2]; line[i+1][1] line[i+1][2]]; digits=digits)
     for i in 1:length(line)-1]
end



function check_hw_fw_all(faults, block_df; verbose=false)

    block_segs = Dict(block.fid => line_to_segs(block.geometry.coords)
                      for block in eachrow(block_df))

    filter(f -> check_hw_fw(f, block_segs; verbose=verbose), faults)
end


function check_hw_fw(fault::Fault, block_segs::Dict; verbose=false)
    
    if !(haskey(block_segs, fault.hw)) || !(haskey(block_segs, fault.fw))
        test_pass = false
    else
        fault_trace = line_to_segs(fault.trace)
        test_pass = check_hw_fw(fault_trace, block_segs[fault.hw], block_segs[fault.fw])
    end

    if !test_pass && verbose
        fid = fault.fid
        @info "removing $fid: problem with hw and/or fw"
    end
    test_pass
end


function check_hw_fw(fault::Fault, hw_poly::Oiler.Geom.Polygon,
    fw_poly::Oiler.Geom.Polygon)

    fault_trace = line_to_segs(fault.trace)
    hw_segs = line_to_segs(hw_poly.coords)
    fw_segs = line_to_segs(fw_poly.coords)

    check_hw_fw(fault_trace, hw_segs, fw_segs)
end


function check_hw_fw(fault_trace, hw_segs, fw_segs)

    hw_adjacent = false
    fw_adjacent = false

    for trace_seg in fault_trace
        if (trace_seg in hw_segs)
            hw_adjacent = true
        end
        if (reverse(trace_seg, dims=1) in fw_segs)
            fw_adjacent = true
        end
    end

    return hw_adjacent && fw_adjacent
end


function get_shared_boundaries(block_df)
    block_seg_sets = Dict(
        block.fid => Set(Oiler.Utils.line_to_segs(block.geometry.coords))
        for block in eachrow(block_df))

    block_seg_sets_rev = Dict(block.fid => Set(Oiler.Utils.line_to_segs(
        reverse(block.geometry.coords, dims=1)))
                              for block in eachrow(block_df))

    block_ids = block_df[!, :fid]

    shared_boundary_segs = Dict()
    for (i, id) in enumerate(block_ids)
        for fid in block_ids[i+1:end]
            shared_segs = intersect(block_seg_sets[id], block_seg_sets_rev[fid])
            if length(shared_segs) > 0
                shared_boundary_segs[(id, fid)] = shared_segs
            end
        end
    end

    shared_boundary_segs
end


function remove_fault_segs_from_bounds(shared_boundary_segs, faults)
    # duplicate faults frwd and rev, don't worry about winding order rn
    shared_boundary_segs_trimmed = Dict()

    for (block_ids, ssegs) in shared_boundary_segs
        shared_boundary_segs_trimmed[block_ids] = ssegs
        for fault in faults
            if (fault.fw in block_ids) & (fault.hw in block_ids)
                fault_segs = Set(Oiler.Utils.line_to_segs(fault.trace))
                fault_segs_r = Set(Oiler.Utils.line_to_segs(
                    reverse(fault.trace, dims=1)))

                ssegs_trimmed = setdiff(shared_boundary_segs_trimmed[block_ids], fault_segs)
                ssegs_trimmed = setdiff(ssegs_trimmed, fault_segs_r)

                shared_boundary_segs_trimmed[block_ids] = ssegs_trimmed

                if length(ssegs) == length(ssegs_trimmed)
                    println("problems with ", block_ids)

                end
            end
        end
        if length(shared_boundary_segs_trimmed[block_ids]) == 0
            delete!(shared_boundary_segs_trimmed, block_ids)
        end
    end
    shared_boundary_segs_trimmed
end


function get_pts_from_segs(segs)
    arr = collect(segs)

    pts = []
    for row in eachrow(arr)
        push!(pts, row[1][1, :])
        push!(pts, row[1][2, :])
    end
    pts
end


function get_endpoints(segs)
    pts = get_pts_from_segs(segs)

    counts = Dict()

    for pt in pts
        if haskey(counts, Tuple(pt))
            counts[Tuple(pt)] += 1
        else
            counts[Tuple(pt)] = 1
        end
    end

    endpts = []
    for (pt, count) in counts
        if count == 1
            push!(endpts, pt)
        end
    end

    endpts

end


function sort_segs(segs; debug=false, max_ends=18, max_segs=1_000)

    segz = collect(segs)

    if debug
        println(length(segs), " segments")
    end

    if length(segz) > max_segs
        return []
    end

    pts = Tuple.(unique(get_pts_from_segs(segs)))

    hash_pts = Dict(pt => i for (i, pt) in enumerate(pts))

    endpts = get_endpoints(segs)

    ind_ends = [hash_pts[e] for e in endpts]
    if debug
        println(length(ind_ends), " endpoints")
    end

    if (length(ind_ends) < 2) || (length(ind_ends) % 2 == 1) || (length(ind_ends) > max_ends)
        return []
    end

    seg_inds = Dict(i => Dict("start" => hash_pts[Tuple(seg[1, :])],
        "end" => hash_pts[Tuple(seg[2, :])])
                    for (i, seg) in enumerate(segz))

    start_seg_inds = Dict(val["start"] => key for (key, val) in seg_inds)

    out_segs = []

    while length(ind_ends) > 0

        first_seg = -1
        for st_ind in ind_ends
            for (si, seg) in seg_inds
                if seg["start"] == st_ind
                    first_seg = si
                    @goto get_going
                end
            end
        end

        @label get_going

        inds = [first_seg]
        this_seg = first_seg

        while !(seg_inds[this_seg]["end"] in ind_ends)
            next_seg = start_seg_inds[seg_inds[this_seg]["end"]]
            push!(inds, next_seg)
            this_seg = next_seg
        end

        coords = zeros(length(inds) + 1, 2)

        for (i, ind) in enumerate(inds)
            coord = pts[seg_inds[ind]["start"]]
            coords[i, 1] = coord[1]
            coords[i, 2] = coord[2]
        end

        coords[end, 1] = pts[seg_inds[last(inds)]["end"]][1]
        coords[end, 2] = pts[seg_inds[last(inds)]["end"]][2]

        push!(out_segs, coords)

        used_ends = [seg_inds[first_seg]["start"], seg_inds[this_seg]["end"]]
        ind_ends = [i for i in ind_ends if !(i in used_ends)]

    end
    out_segs
end

function sort_bound_sets(bound_sets; debug=false)
    dd = Dict()

    n_bound_pairs = length(keys(bound_sets))

    if debug
        println("$n_bound_pairs boundary sets")
    end


    for (i, (pair, segs)) in enumerate(bound_sets)
        try
            if debug
                println("$i / $n_bound_pairs : $pair")
            end
            dd[pair] = Oiler.Utils.sort_segs(segs; debug=debug)
        catch e
            warn_msg = "Can't sort $pair"
            @warn warn_msg
        end
    end
    dd
end


end # module

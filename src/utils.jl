module Utils

export predict_vels_from_poles, find_vel_cycles, diagonalize_matrices, 
    get_gnss_vels, get_coords_from_vel_array

using ..Oiler
using ..Oiler: VelocityVectorSphere, PoleCart, PoleSphere, build_PvGb_from_vels,
        build_vel_column_from_vels, add_poles, pole_sphere_to_cart

using Logging
using DataFrames
using SparseArrays
using LinearAlgebra


function diagonalize_matrices(matrices)

    rowz = [size(m, 1) for m in matrices]
    colz = [size(m, 2) for m in matrices]

    n_rows = sum(rowz)
    n_cols = sum(colz)

    big_mat = spzeros(n_rows, n_cols)

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
function lin_indep_cols(X; tol = 1e-10)
    if ~(true) # supposed to check for all zeros here
        idx = [];
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


function lin_indep_rows(X; tol = 1e-10)
    lin_indep_cols(X'; tol = tol)
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
            @warn "$vel has same fix and mov, leaving it out."
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


function find_tricycles(graph::Dict; reduce::Bool = true)

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

    v_ind = findall(x->x == cycle_tup, vel_group_keys)
    if v_ind != []
        res = Dict("ind" => v_ind[1], "val" => 1.)
    else
        cycle_tup = Base.reverse(cycle_tup)
        v_ind = findall(x->x == (cycle_tup), vel_group_keys)
        if v_ind != []
            res = Dict("ind" => v_ind[1], "val" => -1.)
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
        for ctup in [(tri[1], tri[2]),(tri[2], tri[3]),(tri[3], tri[1])]
            cycle_inds[i][ctup] = get_cycle_inds(vel_group_keys, ctup)
        end
    end
    cycle_inds
end


"""
    flat(x)
Flattens nested arrays of arrays
"""
flat(x, y = vcat(x...)) = x == y ? x : flat(y)


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

    for i in 1:steps
        place = path[i]
        next = path[i + 1]

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
        final_pole = PoleCart(x = 0., y = 0., z = 0., fix = fix, mov = mov)

    elseif length([p for p in poles if (p.fix == fix) & (p.mov == mov)]) == 1
        final_pole = [p for p in poles if (p.fix == fix) & (p.mov == mov)][1]
    elseif length([p for p in poles if (p.mov == fix) & (p.fix == mov)]) == 1
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


function get_vel_vec_at_pole(vel::VelocityVectorSphere, pole::PoleCart)
    PvGb = Oiler.BlockRotations.build_PvGb_vel(vel)
    vel_vec = PvGb * [pole.x; pole.y; pole.z]
end


function get_vel_vec_at_pole(vel::VelocityVectorSphere, pole::PoleSphere)
    get_vel_vec_at_pole(vel, Oiler.BlockRotations.pole_sphere_to_cart(pole))
end


function check_vel_closures(poles; tol = 1e-5)
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
            warn_msg = "$cycle does not close"
            all_good = false
            @warn warn_msg
        end
    end

    if all_good == true
        @info "All poles close"
    end
    all_good
end


function get_fault_slip_rate_from_pole(fault, pole)
    fv = Oiler.Faults.fault_to_vel(fault)
    vel_vec = get_vel_vec_at_pole(fv, pole)
    R = Oiler.Faults.build_strike_rot_matrix(fault.strike)
    v_rl, v_ex, v_up = R * vel_vec
    v_rl, v_ex
end


function get_fault_slip_rates_from_poles(faults, poles)
    rates = []
    for fault in faults
        fault_key = (fault.fw, fault.hw)
        if haskey(poles, fault_key)
            pole = poles[fault_key]
            push!(rates, get_fault_slip_rate_from_pole(fault, pole))
        elseif haskey(poles, reverse(fault_key))
            pole = poles[reverse(fault_key)]
            push!(rates, get_fault_slip_rate_from_pole(fault, pole))
        else
            @warn "No pole found for $fault_key"
            push!(rates, (NaN, NaN))
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
                vel_idx = [row_idx:row_idx + 2, col_idx:col_idx + 2]
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
                vel_idx = [row_idx:row_idx + 2, col_idx:col_idx + 2]
                vd = Dict()
                vd["vel"] = vel
                vd["idx"] = vel_idx
                push!(gnss_vels, vd)
            end # if
        end # for
    end # for
    gnss_vels
end


"""
    get_coords_from_vel_array(vels)

Returns (lons, lats) as arrays of floats from an array of `VelocityVectorSphere`
"""
function get_coords_from_vel_array(vels::Array{VelocityVectorSphere})
    lats = [v.lat for v in vels]
    lons = [v.lon for v in vels]

    (lons, lats)
end

end # module

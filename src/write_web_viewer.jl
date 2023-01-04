module WebViewer

import Base.Threads.@threads

using Mustache
#using ArchGDAL

using ..Oiler

function write_web_viewer(; results, block_df, directory::AbstractString, ref_pole,
    block_filename::AbstractString="blocks.geojson",
    pole_filename::AbstractString="poles.csv",
    c_lon::Float64=110.0, c_lat::Float64=45.0,
    use_center::Bool=true,
    time_min=-5.0, time_max=5.0, time_step=0.25)

    template_file = joinpath(@__DIR__, "web_viewer_template", "web_viewer_template.html")
    blockrot_file = joinpath(@__DIR__, "web_viewer_template", "js", "blockrotations.js")
    main_template_file = joinpath(@__DIR__, "web_viewer_template", "js", "main_template.js")
    main_file = joinpath(@__DIR__, "web_viewer_template", "js", "main.js")

    index_html_template = Mustache.load(template_file)
    main_js_template = Mustache.load(main_template_file)

    block_filepath = joinpath(directory, block_filename)
    pole_filepath = joinpath(directory, pole_filename)

    filled_index_html = index_html_template(time_min=time_min, time_max=time_max,
        time_step=time_step,
        blocks_path=block_filename,
        poles_path=pole_filename)

    get_block_colors(block_df)

    if use_center
        block_centroids = Oiler.Utils.get_block_centroids(block_df)
        c_lon = Oiler.Geom.angular_mean_degrees([c.coords[1] for c in block_centroids])
        c_lat = Oiler.Geom.angular_mean_degrees([c.coords[2] for c in block_centroids])
    end

    main_js_centered = main_js_template(center_lon=-c_lon, center_lat=-c_lat)

    mkpath(joinpath(directory, "js"))
    cp(blockrot_file, joinpath(directory, "js", "blockrotations.js"); force=true)
    #cp(main_file, "joinpath(directory, "js", "main.js")"; force = true)

    Oiler.IO.write_solution_poles(pole_filepath, results, block_df, ref_pole)
    Oiler.IO.write_block_df(block_df, block_filepath, simplify=true)

    open(joinpath(directory, "index.html"), "w") do f
        write(f, filled_index_html)
    end
    open(joinpath(directory, "js", "main.js"), "w") do f
        write(f, main_js_centered)
    end
end


function get_block_adj(block_df)
    df_adj = Dict(block_df[i, :fid] => [] for i in 1:size(block_df, 1))

    @threads for i in 1:size(block_df, 1)
        f1 = block_df[i, :]
        fid1 = f1[:fid]
        g1 = f1[:geometry]
        for j in 1:size(block_df, 1)
            f2 = block_df[j, :]
            fid2 = f2[:fid]
            were_touching = check_geom(g1, f2[:geometry])
            if were_touching
                if !(fid1 in df_adj[fid2])
                    push!(df_adj[fid2], fid1)
                end
                if !(fid2 in df_adj[fid1])
                    push!(df_adj[fid1], fid1)
                end
            end
        end
    end
    df_adj
end


function check_geom(g1, g2; n_min_pts=2)
    g1c = g1.coords
    g2c = g2.coords

    pt_set_1 = Set([Tuple(g1c[i, :]) for i in 1:size(g1c, 1)])
    pt_set_2 = Set([Tuple(g2c[i, :]) for i in 1:size(g2c, 1)])

    if pt_set_1 == pt_set_2
        return_val = false
    else
        common_pts = intersect(pt_set_1, pt_set_2)
        return_val = (length(common_pts) >= n_min_pts)
    end
    return_val
end



const colors = ["#c994c7", "#ccebc5", "#ff5112", "#7bccc4", "#4eb3d3",
    "#2b8cbe", "#08589e", "#fde0dd", "#fcc5c0", "#fa9fb5",
    "#f768a1", "#dd3497", "#ae017e", "#7a0177"]


function get_block_colors(block_df)

    idx_fid = Dict(i => block_df[i, :fid] for i in 1:size(block_df, 1))
    fid_idx = Dict(fid => i for (i, fid) in idx_fid)

    block_colors = ["" for i in 1:size(block_df, 1)]
    block_adj = get_block_adj(block_df)

    for i in 1:size(block_df, 1)
        adj_colors = []
        f = block_df[i, :]
        for adj_fid in block_adj[idx_fid[i]]
            ac = block_colors[fid_idx[adj_fid]]
            if ac != ""
                push!(adj_colors, ac)
            end
        end
        avail_colors = collect(setdiff(Set(colors), Set(adj_colors)))
        if length(avail_colors) == 0
            block_colors[i] = rand(colors)
        else
            block_colors[i] = rand(avail_colors)
        end
    end
    block_df[!, :color] = block_colors
end




end # module

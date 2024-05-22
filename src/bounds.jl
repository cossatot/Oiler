module Boundaries

using ..Oiler

struct Boundary
    trace::Array{Float64,2}
    fix::String
    mov::String
    fid::String
end



function boundary_to_vels(boundary::Boundary; ve=0.0, vn=0.0, ee=1.0, en=1.0,
    simp_dist=2.0)

    simp_trace = Oiler.Geom.simplify_polyline(boundary.trace, simp_dist)
    trace_length = Oiler.Geom.polyline_length(simp_trace)

    if trace_length < 20.0
        n_segs = 1
    elseif (20.0 <= trace_length) & (trace_length < 60.0)
        n_segs = 2
    elseif (60.0 <= trace_length) & (trace_length < 120.0)
        n_segs = 3
    else
        n_segs = Int(floor(trace_length / 30.0))
    end

    new_traces = Oiler.Geom.break_polyline_equal(simp_trace, n_segs)

    midpoints = map(Oiler.Faults.get_midpoint, new_traces)

    vels = map(midpt -> Oiler.VelocityVectorSphere(
            lon=midpt[1],
            lat=midpt[2],
            ve=ve,
            vn=vn,
            ee=ee,
            en=en,
            fix=boundary.fix,
            mov=boundary.mov,
            vel_type="Boundary",
        ),
        midpoints
    )

    vels
end


end #module

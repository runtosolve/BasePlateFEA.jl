# --------------------------------------------------------------------------------------------
# 3D viewer: a self-contained HTML page (three.js) with the deformed model, element stress and
# plastic strain contours, the contact stress on the plate, the anchor bolt forces, and the
# moment–rotation and axial force–shortening curves, for every stored load step.
# --------------------------------------------------------------------------------------------

const VIEWER_TEMPLATE = joinpath(@__DIR__, "..", "viewer", "template.html")

_r(x::Real, n) = isfinite(x) ? round(x; sigdigits = n) : nothing
_r(v::AbstractVector, n) = [_r(x, n) for x in v]

"""
    viewer_data(result; title, subtitle = "", curves = [], max_frames = 40, fields)

Data of the 3D viewer as a `Dict`. `curves` are extra moment–rotation curves drawn with the model's:
named tuples `(label, kind, θ, M)` with `kind = :test` or `:reference`. At most `max_frames` load
steps are exported (evenly thinned; the unloaded state, the end of the axial step and the last
step are always kept). `second_order` reports the base moment including `N δ`
(see [`base_moment`](@ref)). `fields` selects the element results (see [`element_fields`](@ref)).
"""
function viewer_data(res::Result; title::AbstractString = "Base plate model", subtitle::AbstractString = "", curves = [],
                     max_frames::Int = 40, second_order::Bool = false, fields = (:vm_top, :vm_bot, :vm_mid, :peeq, :sx, :sy, :sx_top, :sy_top))
    m = res.model.mesh
    nf = length(res.frames)
    i1 = findlast(f -> f.step == 1, res.frames)
    keep = sort!(unique(vcat(1, i1, nf, round.(Int, range(1, nf; length = min(nf, max_frames))))))
    frames = res.frames[keep]

    tris = Int[]; tri_elem = Int[]; tri_part = Int[]; es = Set{NTuple{2,Int}}()
    for (e, c) in enumerate(m.conn)
        m.hole[e] && continue
        for k in (length(c) == 3 ? ((1, 2, 3),) : ((1, 2, 3), (1, 3, 4)))
            append!(tris, (c[k[1]] - 1, c[k[2]] - 1, c[k[3]] - 1)); push!(tri_elem, e - 1); push!(tri_part, m.part[e])
        end
        for k in eachindex(c)
            a = c[k]; b = c[mod1(k + 1, length(c))]
            push!(es, a < b ? (a - 1, b - 1) : (b - 1, a - 1))
        end
    end
    lo = vec(minimum(m.x, dims = 1)); hi = vec(maximum(m.x, dims = 1))
    size = maximum(hi - lo)
    nn = nnodes(m)
    Mb(f) = base_moment(res.model, f; second_order)
    disp(f) = vec(reshape(f.u[1:6nn], 6, nn)[1:3, :])
    umax = maximum(maximum(abs, disp(f)) for f in frames)
    default_scale = umax > 0 ? clamp(round(0.04 * size / umax; sigdigits = 1), 1.0, 500.0) : 1.0

    fr = [Dict(:step => f.step, :N => f.N, :M => Mb(f), :theta => f.θ, :w => f.w, :bolt => _r(f.bolt_force, 5),
               :contact => _r(f.contact_stress, 4), :u => _r(disp(f), 5),
               :fields => Dict(k => _r(f.fields[k], 4) for k in fields)) for f in frames]
    i1k = findlast(g -> g.step == 1, frames)
    s2 = collect(i1k:length(frames))
    return Dict(
        :title => title, :subtitle => subtitle,
        :nodes => _r(vec(permutedims(m.x)), 7), :tris => tris, :tri_elem => tri_elem, :tri_part => tri_part,
        :edges => reduce(vcat, [[a, b] for (a, b) in es]; init = Int[]),
        :bp_nodes => m.bp_nodes .- 1, :bolt_nodes => m.bolt_nodes .- 1,
        :center => (lo + hi) / 2, :size => size, :default_scale => default_scale,
        :frames => fr,
        :mr => Dict(:theta => [frames[i].θ for i in s2], :M => [Mb(frames[i]) for i in s2], :frame => s2 .- 1),
        :nd => Dict(:w => [-f.w for f in frames], :N => [f.N for f in frames], :frame => collect(0:length(frames)-1)),
        :curves => [Dict(:label => c.label, :kind => String(c.kind), :x => _r(c.θ, 5), :y => _r(c.M, 5)) for c in curves],
    )
end

"""
    write_viewer(path, result; kwargs...)

Write the 3D viewer of `result` to the HTML file `path` (keywords as for [`viewer_data`](@ref)).
The page loads three.js from a CDN, so it needs an internet connection when opened.
"""
function write_viewer(path::AbstractString, res::Result; kwargs...)
    html = read(VIEWER_TEMPLATE, String)
    occursin("/*DATA*/null", html) || error("viewer template has no data placeholder")
    data = replace(JSON3.write(viewer_data(res; kwargs...)), "</" => "<\\/")
    write(path, replace(html, "/*DATA*/null" => data))
    return path
end

# --------------------------------------------------------------------------------------------
# Makie / Bonito 3D viewer of BasePlateFEA results.
#
#   serve_viewer(result; ...)          live page served by Bonito (WGLMakie)
#   export_viewer(path, result; ...)   standalone HTML page (no Julia needed to view it)
#   render_image(path, result, frame)  high-quality still image (GLMakie: SSAO, FXAA, supersampling)
#
# The shells are drawn with their thickness (two offset surfaces and the free edges closed), with
# banded contours of nodal values averaged within each part (base plate, each upright wall).
# All per-step data are precomputed and sent to the browser once; the step slider, field selector,
# deformation scale, contour limit and toggles are handled in JavaScript, so the live page and the
# exported page are the same application.
# --------------------------------------------------------------------------------------------
module BasePlateViewer

using BasePlateFEA
using Makie, Bonito, GeometryBasics, Colors, Printf, LinearAlgebra
import WGLMakie, GLMakie

export serve_viewer, export_viewer, render_image, viewer_app, VizData

const NBANDS = 12
const QNAN = typemax(UInt16)

# ---- palettes ------------------------------------------------------------------------------

_hex(s) = parse(RGBf, s)
function _ramp(stops::Vector{RGBf}, n)
    out = RGBf[]
    for i in 1:n
        x = (i - 0.5) / n * (length(stops) - 1)
        k = clamp(floor(Int, x), 0, length(stops) - 2)
        f = Float32(x - k)
        a = stops[k+1]; b = stops[k+2]
        push!(out, RGBf(a.r + f * (b.r - a.r), a.g + f * (b.g - a.g), a.b + f * (b.b - a.b)))
    end
    return out
end

"`NBANDS` colors for magnitudes: `:abaqus` (blue–green–red spectrum, default), `:viridis`, `:blue` (single hue)"
function sequential_bands(palette::Symbol)
    palette == :abaqus && return [RGBf(HSV(240 * (1 - (i - 1) / (NBANDS - 1)), 0.82, 0.96)) for i in 1:NBANDS]
    palette == :blue && return _ramp(_hex.(["#e4effd", "#9ec5f4", "#3987e5", "#184f95", "#0a2a55"]), NBANDS)
    palette == :viridis && return [RGBf(c) for c in Makie.to_colormap(cgrad(:viridis, NBANDS; categorical = true))]
    throw(ArgumentError("palette must be :viridis, :abaqus or :blue"))
end
"`NBANDS` colors for signed values (blue – neutral – orange); `:abaqus` keeps the spectrum"
diverging_bands(palette::Symbol) = palette == :abaqus ? sequential_bands(:abaqus) :
                                   _ramp(_hex.(["#0b3168", "#3987e5", "#f0efec", "#eb6834", "#7c2c0a"]), NBANDS)

# ---- precomputed data ----------------------------------------------------------------------

struct FieldDef
    key::Symbol
    label::String
    unit::String
    kind::Symbol          # :seq (0 … max) or :div (−max … max)
    top::Symbol           # stored array shown on the +normal face
    bot::Symbol           # stored array shown on the −normal face
    auto::Float64         # default contour limit
    plate_only::Bool
end

"""
    VizData(result; max_frames = 40, second_order = false, thickness_scale = 1.0)

Everything the viewer draws, precomputed from a `BasePlateFEA.Result`.
"""
struct VizData
    x::Matrix{Float32}                 # nn × 3 node coordinates
    pv_node::Vector{Int32}             # value slots: one per (part, node)
    pv_part::Vector{Int32}
    pv_off::Matrix{Float32}            # npv × 3 half-thickness offset vector (+normal)
    sv_pv::Vector{Int32}               # surface vertices → value slot
    sv_side::Vector{Int8}              # +1 top (+normal) face, −1 bottom face
    faces::Vector{NTuple{3,Int32}}     # triangles of surface vertices (1-based)
    face_part::Vector{Int32}
    edges::Vector{NTuple{2,Int32}}     # element edges (surface vertices), both faces
    edge_part::Vector{Int32}
    feature::Vector{NTuple{2,Int32}}   # part boundaries and thickness edges
    feature_part::Vector{Int32}
    bolt_nodes::Vector{Int32}
    disp::Array{Float32,3}             # 3 × nn × nframes
    stored::Dict{Symbol,Matrix{UInt16}}       # npv × nframes, quantized; QNAN = no value
    range::Dict{Symbol,NTuple{2,Float64}}     # value range of each stored array
    fields::Vector{FieldDef}
    step::Vector{Int}; N::Vector{Float64}; M::Vector{Float64}; θ::Vector{Float64}; w::Vector{Float64}
    bolt::Matrix{Float64}              # nbolts × nframes
    contact_max::Vector{Float64}
    first2::Int                        # index of the last frame of the axial step
    center::Vector{Float64}; size::Float64; default_scale::Float64
    plate_t::Float64
end

function _nice(v)
    (isfinite(v) && v > 0) || return 1.0
    e = 10.0^floor(log10(v)); f = v / e
    return (f <= 1 ? 1.0 : f <= 1.2 ? 1.2 : f <= 1.5 ? 1.5 : f <= 2 ? 2.0 : f <= 2.5 ? 2.5 : f <= 3 ? 3.0 : f <= 4 ? 4.0 :
            f <= 5 ? 5.0 : f <= 6 ? 6.0 : f <= 8 ? 8.0 : 10.0) * e
end

function VizData(res::BasePlateFEA.Result; max_frames::Int = 40, second_order::Bool = false, thickness_scale::Real = 1.0)
    m = res.model.mesh
    nn = size(m.x, 1)
    nf_all = length(res.frames)
    i1 = findlast(f -> f.step == 1, res.frames)
    keep = sort!(unique(vcat(1, i1, nf_all, round.(Int, range(1, nf_all; length = min(nf_all, max_frames))))))
    frames = res.frames[keep]
    nf = length(frames)
    els = [e for e in 1:length(m.conn) if !m.hole[e]]

    # value slots: one per (part, node)
    slot = Dict{Tuple{Int,Int},Int32}(); pv_node = Int32[]; pv_part = Int32[]; pv_t = Float64[]
    for e in els, i in m.conn[e]
        get!(slot, (m.part[e], i)) do
            push!(pv_node, i); push!(pv_part, m.part[e]); push!(pv_t, m.thickness[e])
            Int32(length(pv_node))
        end
    end
    npv = length(pv_node)
    pv_off = zeros(Float32, npv, 3)
    for v in 1:npv
        pv_off[v, :] .= thickness_scale * pv_t[v] / 2 .* m.frames[pv_part[v]][:, 3]
    end
    nparts_of = zeros(Int, nn)
    for v in 1:npv
        nparts_of[pv_node[v]] += 1
    end

    # surface vertices: top copy (1…npv), bottom copy (npv+1…2npv), then the closing strips
    sv_pv = Int32[1:npv; 1:npv]; sv_side = Int8[fill(Int8(1), npv); fill(Int8(-1), npv)]
    faces = NTuple{3,Int32}[]; face_part = Int32[]
    ecount = Dict{NTuple{2,Int32},Int}(); epart = Dict{NTuple{2,Int32},Int32}()
    for e in els
        p = m.part[e]
        c = [slot[(p, i)] for i in m.conn[e]]
        for k in (length(c) == 3 ? ((1, 2, 3),) : ((1, 2, 3), (1, 3, 4)))
            push!(faces, (c[k[1]], c[k[2]], c[k[3]])); push!(face_part, p)
            push!(faces, (c[k[1]] + npv, c[k[3]] + npv, c[k[2]] + npv)); push!(face_part, p)
        end
        for k in eachindex(c)
            a = c[k]; b = c[mod1(k + 1, length(c))]
            key = a < b ? (a, b) : (b, a)
            ecount[key] = get(ecount, key, 0) + 1; epart[key] = p
        end
    end
    edges = NTuple{2,Int32}[]; edge_part = Int32[]; feature = NTuple{2,Int32}[]; feature_part = Int32[]
    for (key, cnt) in ecount
        a, b = key; p = epart[key]
        push!(edges, (a, b)); push!(edges, (a + npv, b + npv)); push!(edge_part, p); push!(edge_part, p)
        cnt == 1 || continue
        push!(feature, (a, b)); push!(feature, (a + npv, b + npv)); push!(feature_part, p); push!(feature_part, p)
        (nparts_of[pv_node[a]] > 1 && nparts_of[pv_node[b]] > 1) && continue     # junction of two parts: not a free edge
        n0 = Int32(length(sv_pv))                                              # free edge: close the thickness
        append!(sv_pv, (a, b, b, a)); append!(sv_side, (1, 1, -1, -1))
        push!(faces, (n0 + 1, n0 + 2, n0 + 3)); push!(faces, (n0 + 1, n0 + 3, n0 + 4)); push!(face_part, p); push!(face_part, p)
    end

    # displacements and stored fields
    disp = zeros(Float32, 3, nn, nf)
    for (k, f) in enumerate(frames)
        disp[:, :, k] .= reshape(f.u[1:6nn], 6, nn)[1:3, :]
    end
    ekeys = (:vm_top, :vm_bot, :vm_mid, :peeq, :sx, :sy, :sx_top, :sx_bot, :sy_top, :sy_bot)
    raw = Dict{Symbol,Matrix{Float64}}(k => fill(NaN, npv, nf) for k in (ekeys..., :contact, :uz))
    bpidx = Dict(n => i for (i, n) in enumerate(m.bp_nodes))
    acc = zeros(npv); cnt = zeros(Int, npv)
    for (k, f) in enumerate(frames)
        for key in ekeys
            v = f.fields[key]
            fill!(acc, 0.0); fill!(cnt, 0)
            for e in els
                isnan(v[e]) && continue
                for i in m.conn[e]
                    s = slot[(m.part[e], i)]
                    acc[s] += v[e]; cnt[s] += 1
                end
            end
            for s in 1:npv
                cnt[s] > 0 && (raw[key][s, k] = acc[s] / cnt[s])
            end
        end
        for s in 1:npv
            raw[:uz][s, k] = f.u[6*(pv_node[s]-1)+3]
            pv_part[s] == 1 && haskey(bpidx, pv_node[s]) && (raw[:contact][s, k] = f.contact_stress[bpidx[pv_node[s]]])
        end
    end
    stored = Dict{Symbol,Matrix{UInt16}}(); rng = Dict{Symbol,NTuple{2,Float64}}()
    for (key, A) in raw
        fin = filter(isfinite, A)
        lo, hi = isempty(fin) ? (0.0, 1.0) : extrema(fin)
        hi > lo || (hi = lo + 1.0)
        rng[key] = (lo, hi)
        stored[key] = map(v -> isfinite(v) ? round(UInt16, (v - lo) / (hi - lo) * 65534) : QNAN, A)
    end
    amax(keys...) = _nice(maximum(max(abs(rng[k][1]), abs(rng[k][2])) for k in keys))
    fields = [
        FieldDef(:vm_surf, "von Mises stress, shell surfaces", "ksi", :seq, :vm_top, :vm_bot, amax(:vm_top, :vm_bot), false),
        FieldDef(:vm_mid, "von Mises stress, mid-surface", "ksi", :seq, :vm_mid, :vm_mid, amax(:vm_mid), false),
        FieldDef(:vm_top, "von Mises stress, top (+normal) surface", "ksi", :seq, :vm_top, :vm_top, amax(:vm_top), false),
        FieldDef(:vm_bot, "von Mises stress, bottom (−normal) surface", "ksi", :seq, :vm_bot, :vm_bot, amax(:vm_bot), false),
        FieldDef(:peeq, "Equivalent plastic strain", "", :seq, :peeq, :peeq, amax(:peeq), false),
        FieldDef(:contact, "Contact stress on the concrete", "ksi", :seq, :contact, :contact, amax(:contact), true),
        FieldDef(:sx_surf, "σx, shell surfaces (plate x / along wall)", "ksi", :div, :sx_top, :sx_bot, amax(:sx_top, :sx_bot), false),
        FieldDef(:sy_surf, "σy, shell surfaces (plate y / vertical)", "ksi", :div, :sy_top, :sy_bot, amax(:sy_top, :sy_bot), false),
        FieldDef(:sx, "σx, membrane", "ksi", :div, :sx, :sx, amax(:sx), false),
        FieldDef(:sy, "σy, membrane", "ksi", :div, :sy, :sy, amax(:sy), false),
        FieldDef(:uz, "Vertical displacement", "in", :div, :uz, :uz, amax(:uz), false),
    ]
    lo = vec(minimum(m.x, dims = 1)); hi = vec(maximum(m.x, dims = 1))
    size_ = maximum(hi - lo)
    umax = maximum(abs, disp)
    default_scale = umax > 0 ? clamp(round(0.04 * size_ / umax; sigdigits = 1), 1.0, 500.0) : 1.0
    i1k = findlast(f -> f.step == 1, frames)
    return VizData(Float32.(m.x), pv_node, pv_part, pv_off, sv_pv, sv_side, faces, face_part, edges, edge_part, feature, feature_part,
                   Int32.(m.bolt_nodes), disp, stored, rng, fields,
                   [f.step for f in frames], [f.N for f in frames], [BasePlateFEA.base_moment(res.model, f; second_order) for f in frames],
                   [f.θ for f in frames], [f.w for f in frames], reduce(hcat, [f.bolt_force for f in frames]),
                   [maximum(f.contact_stress) for f in frames], i1k, (lo + hi) / 2, size_, default_scale, res.model.inp.geometry.t)
end

nframes(v::VizData) = length(v.step)
fieldindex(v::VizData, key::Symbol) = something(findfirst(f -> f.key == key, v.fields))

# ---- per-frame arrays (the JavaScript in `viewer_app` mirrors these) --------------------------

function surface_positions(v::VizData, k::Int, scale::Real)
    out = Vector{Point3f}(undef, length(v.sv_pv))
    for i in eachindex(v.sv_pv)
        s = v.sv_pv[i]; n = v.pv_node[s]; sg = v.sv_side[i]
        out[i] = Point3f(v.x[n, 1] + scale * v.disp[1, n, k] + sg * v.pv_off[s, 1],
                         v.x[n, 2] + scale * v.disp[2, n, k] + sg * v.pv_off[s, 2],
                         v.x[n, 3] + scale * v.disp[3, n, k] + sg * v.pv_off[s, 3])
    end
    return out
end

# band coordinate in [0, 2 NBANDS): sequential fields use the first NBANDS colors, diverging the last
function surface_colors(v::VizData, k::Int, fd::FieldDef, cmax::Real)
    out = Vector{Float32}(undef, length(v.sv_pv))
    lo, hi = fd.kind == :seq ? (0.0, float(cmax)) : (-float(cmax), float(cmax))
    base = fd.kind == :seq ? 0.0 : float(NBANDS)
    for i in eachindex(v.sv_pv)
        key = v.sv_side[i] > 0 ? fd.top : fd.bot
        q = v.stored[key][v.sv_pv[i], k]
        if q == QNAN
            out[i] = NaN32
        else
            r = v.range[key]
            val = r[1] + q / 65534 * (r[2] - r[1])
            out[i] = base + NBANDS * clamp((val - lo) / (hi - lo), 0.0, 1 - 1e-6)
        end
    end
    return out
end

segment_points(pos, segs) = Point3f[pos[i] for s in segs for i in s]
function bolt_positions(v::VizData, k::Int, scale::Real)
    return [Point3f(v.x[n, 1] + scale * v.disp[1, n, k], v.x[n, 2] + scale * v.disp[2, n, k],
                    v.x[n, 3] + scale * v.disp[3, n, k] + v.plate_t / 2) for n in v.bolt_nodes]
end

band_colormap(palette) = cgrad(vcat(sequential_bands(palette), diverging_bands(palette)), 2NBANDS; categorical = true)

# ---- scene ---------------------------------------------------------------------------------

"draw the model for frame `k` into the 3D scene `ax`; returns the plots"
function draw_model!(ax, v::VizData, k::Int, fd::FieldDef; scale, cmax, palette = :abaqus, edges = true, upright = true,
                     edge_width = 0.6, feature_width = 1.6, ssao = false)
    pos = surface_positions(v, k, scale)
    fsel = upright ? eachindex(v.faces) : findall(==(1), v.face_part)
    faces = [GLTriangleFace(f...) for f in v.faces[fsel]]
    msh = mesh!(ax, pos, faces; color = surface_colors(v, k, fd, cmax), colormap = band_colormap(palette), colorrange = (0, 2NBANDS),
                nan_color = RGBf(0.80, 0.80, 0.78), shading = true, backlight = 1.0f0, specular = 0.12f0, shininess = 24.0f0,
                diffuse = 0.55f0, interpolate = false, ssao)
    esel = upright ? eachindex(v.edges) : findall(==(1), v.edge_part)
    ed = linesegments!(ax, segment_points(pos, v.edges[esel]); color = RGBAf(0.05, 0.05, 0.05, 0.38), linewidth = edge_width,
                       visible = edges, depth_shift = -1.0f-4)
    gsel = upright ? eachindex(v.feature) : findall(==(1), v.feature_part)
    ft = linesegments!(ax, segment_points(pos, v.feature[gsel]); color = RGBf(0.04, 0.04, 0.04), linewidth = feature_width, depth_shift = -3.0f-4)
    pos0 = surface_positions(v, 1, 0.0)
    gh = linesegments!(ax, segment_points(pos0, v.feature); color = RGBAf(0.35, 0.35, 0.33, 0.6), linewidth = 0.8, visible = false)
    R = 0.3f0
    bl = meshscatter!(ax, bolt_positions(v, k, scale); marker = GeometryBasics.Cylinder(Point3f(0, 0, 0), Point3f(0, 0, 1), 1.0f0),
                      markersize = Vec3f(R, R, 0.18), color = RGBf(0.62, 0.64, 0.67), shading = true, specular = 0.3f0, ssao)
    return (; mesh = msh, edges = ed, feature = ft, ghost = gh, bolts = bl)
end

function set_camera!(lscene, v::VizData; azimuth = -0.95, elevation = 0.42, zoom = 1.0)
    c = Vec3f(v.center...)
    d = Vec3f(cos(elevation) * cos(azimuth), cos(elevation) * sin(azimuth), sin(elevation))
    cam = cameracontrols(lscene.scene)
    cam.fov[] = 28.0f0
    update_cam!(lscene.scene, c + Float32(1.65 * v.size / zoom) * d, c, Vec3f(0, 0, 1))
    return lscene
end

_lights() = [AmbientLight(RGBf(0.80, 0.80, 0.80)), DirectionalLight(RGBf(0.36, 0.36, 0.36), Vec3f(-0.35, 0.45, -1.0), true)]

function charts!(gl, v::VizData, k::Int; curves = [])
    s1 = RGBf(_hex("#2a78d6")); s2 = RGBf(_hex("#eb6834")); ink = RGBf(0.32, 0.32, 0.31); grid = RGBf(0.88, 0.88, 0.85)
    kw = (; xgridcolor = grid, ygridcolor = grid, topspinevisible = false, rightspinevisible = false, titlealign = :left,
          titlesize = 14, titlefont = :regular, xlabelsize = 12, ylabelsize = 12, xticklabelsize = 11, yticklabelsize = 11,
          xlabelcolor = ink, ylabelcolor = ink, xticklabelcolor = ink, yticklabelcolor = ink)
    a1 = Axis(gl[1, 1]; title = "Base moment – rotation", xlabel = "Rotation (rad)", ylabel = "Moment (kip-in)", kw...)
    for c in curves
        String(c.kind) == "test" && lines!(a1, c.θ, c.M; color = RGBAf(0.35, 0.35, 0.33, 0.55), linewidth = 1.0)
    end
    for c in curves
        String(c.kind) == "reference" && lines!(a1, c.θ, c.M; color = s2, linewidth = 1.6, linestyle = :dash, label = c.label)
    end
    r2 = v.first2:nframes(v)
    lines!(a1, v.θ[r2], v.M[r2]; color = s1, linewidth = 2.2, label = "Model")
    any(c -> String(c.kind) == "test", curves) && lines!(a1, [NaN], [NaN]; color = RGBAf(0.35, 0.35, 0.33, 0.8), linewidth = 1.0, label = "Tests")
    m1 = scatter!(a1, [Point2f(v.θ[k], v.M[k])]; color = s1, strokecolor = :white, strokewidth = 2, markersize = 13)
    isempty(curves) || axislegend(a1; position = :rb, framevisible = false, labelsize = 11, patchsize = (18, 10), rowgap = 0, padding = (4, 4, 2, 2))
    a2 = Axis(gl[2, 1]; title = "Axial force – shortening", xlabel = "Shortening at the rotation gauges (in)", ylabel = "Axial force (kip)", kw...)
    r1 = 1:v.first2
    lines!(a2, -v.w[r1], v.N[r1]; color = s1, linewidth = 2.2)
    lines!(a2, -v.w[r2], v.N[r2]; color = s1, linewidth = 1.2, linestyle = :dot)
    m2 = scatter!(a2, [Point2f(-v.w[k], v.N[k])]; color = s1, strokecolor = :white, strokewidth = 2, markersize = 13)
    return m1, m2
end

fmt(x) = abs(x) >= 100 ? @sprintf("%.0f", x) : abs(x) >= 10 ? @sprintf("%.1f", x) : abs(x) >= 1 ? @sprintf("%.2f", x) :
         x == 0 ? "0" : @sprintf("%.3g", x)

# ---- still images (GLMakie) ----------------------------------------------------------------

"""
    render_image(path, result_or_vizdata, frame = last; field = :vm_surf, scale = auto, cmax = auto,
                 palette = :abaqus, size = (1600, 1000), px_per_unit = 2, ssao = true, edges = true,
                 upright = true, azimuth, elevation, zoom, title = "", curves = [], charts = false)

Write a PNG of load step `frame` (index into the exported frames) with GLMakie: banded contours with a
legend, element edges, screen-space ambient occlusion and anti-aliasing at `px_per_unit` supersampling.
"""
function render_image(path::AbstractString, v::VizData, frame::Int = nframes(v); field::Symbol = :vm_surf, scale = v.default_scale,
                      cmax = nothing, palette::Symbol = :abaqus, size = (1600, 1000), px_per_unit = 2, ssao::Bool = true,
                      edges::Bool = true, upright::Bool = true, azimuth = -0.95, elevation = 0.42, zoom = 1.0, title::AbstractString = "",
                      curves = [], charts::Bool = false)
    GLMakie.activate!(; ssao, fxaa = true, px_per_unit, visible = false)
    fd = v.fields[fieldindex(v, field)]
    cm = cmax === nothing ? fd.auto : float(cmax)
    fig = Figure(; size, backgroundcolor = :white, figure_padding = 14)
    ls = LScene(fig[1, 1]; show_axis = false, scenekw = (; lights = _lights(), backgroundcolor = RGBf(0.96, 0.96, 0.95), clear = true))
    draw_model!(ls, v, frame, fd; scale, cmax = cm, palette, edges, upright, ssao)
    set_camera!(ls, v; azimuth, elevation, zoom)
    bands = fd.kind == :seq ? sequential_bands(palette) : diverging_bands(palette)
    lo = fd.kind == :seq ? 0.0 : -cm
    ticks = collect(range(lo, cm; length = NBANDS + 1))
    Colorbar(fig[1, 2]; colormap = cgrad(bands, NBANDS; categorical = true), limits = (lo, cm), ticks = (ticks, fmt.(ticks)),
             label = isempty(fd.unit) ? fd.label : "$(fd.label) ($(fd.unit))", width = 22, ticklabelsize = 13, labelsize = 14,
             height = Relative(0.6), tellheight = false)
    info = @sprintf("N = %.2f kip    M = %.2f kip-in    θ = %.4f rad    deformation ×%s", v.N[frame], v.M[frame], v.θ[frame], fmt(scale))
    Label(fig[0, 1:2], isempty(title) ? info : title * "\n" * info; halign = :left, fontsize = 15, justification = :left)
    if charts
        gl = GridLayout(fig[1, 3]); charts!(gl, v, frame; curves); colsize!(fig.layout, 3, Fixed(430))
    end
    save(path, fig; px_per_unit)
    WGLMakie.activate!()
    return path
end
render_image(path::AbstractString, res::BasePlateFEA.Result, frame::Int = 0; max_frames::Int = 40, second_order::Bool = false, kw...) =
    (v = VizData(res; max_frames, second_order); render_image(path, v, frame == 0 ? nframes(v) : frame; kw...))

# ---- Bonito application ----------------------------------------------------------------------

const CSS = """
body { margin: 0; background: #f9f9f7; color: #0b0b0b; font: 14px/1.4 -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif; }
.bpv { max-width: 1680px; margin: 0 auto; padding: 14px 18px 22px; }
.bpv h1 { font-size: 20px; font-weight: 600; margin: 0; }
.bpv .sub { color: #52514e; margin: 2px 0 10px; }
.bpv .bar { display: flex; flex-wrap: wrap; gap: 8px 18px; align-items: center; padding: 8px 12px; background: #fcfcfb;
            border: 1px solid #e1e0d9; border-radius: 8px; margin-bottom: 10px; }
.bpv .bar label { color: #52514e; font-size: 13px; display: flex; align-items: center; gap: 6px; }
.bpv select, .bpv input[type=number], .bpv button { font: inherit; font-size: 13px; padding: 3px 6px; border: 1px solid #c3c2b7;
            border-radius: 5px; background: #fff; color: #0b0b0b; }
.bpv input[type=number] { width: 72px; }
.bpv button { cursor: pointer; min-width: 58px; }
.bpv .step { flex: 1 1 320px; }
.bpv .step input { flex: 1; }
.bpv .wrap { position: relative; display: flex; align-items: stretch; background: #fff; border: 1px solid #e1e0d9; border-radius: 8px; overflow: hidden; }
.bpv .left { flex: 0 0 168px; padding: 14px 0 0 14px; }
.bpv .right { flex: 1 1 auto; min-width: 0; position: relative; padding-bottom: 30px; }
.bpv .legend { font-size: 12px; font-variant-numeric: tabular-nums; }
.bpv .legend .ttl { font-weight: 600; margin-bottom: 14px; max-width: 150px; }
.bpv .legend .row { display: flex; align-items: flex-start; height: 22px; }
.bpv .legend .sw { width: 30px; height: 22px; box-sizing: border-box; border: 1px solid rgba(0,0,0,0.55); border-top-width: 0; }
.bpv .legend .row:first-of-type .sw { border-top-width: 1px; }
.bpv .legend .lab { margin-left: 7px; transform: translateY(-8px); color: #0b0b0b; }
.bpv .legend .last { height: 0; }
.bpv .legend .note { color: #52514e; margin-top: 12px; }
.bpv .info { position: absolute; left: 8px; bottom: 8px; font-size: 12.5px; font-variant-numeric: tabular-nums; }
.bpv .info b { font-weight: 600; }
"""

_u16(A) = vec(A)
_i32(v) = Int32.(v)

"""
    viewer_app(result; title, subtitle = "", curves = [], max_frames = 40, second_order = false,
               palette = :abaqus, size = (1420, 820), frame = last, field = :vm_surf) -> Bonito.App

The viewer as a Bonito application (see [`serve_viewer`](@ref), [`export_viewer`](@ref)).
`curves` are extra moment–rotation curves, named tuples `(label, kind, θ, M)` with `kind = :test`
or `:reference`. `second_order` reports the base moment including `N δ` (`BasePlateFEA.base_moment`).
`palette` is `:abaqus` (spectrum), `:viridis` or `:blue`; signed fields use a blue–orange diverging
scale except with `:abaqus`. The page accepts an initial state in the URL hash, e.g.
`#frame=12&field=peeq&scale=5&cmax=0.02&upright=0&edges=0`.
"""
function viewer_app(v::VizData; title::AbstractString = "Base plate model", subtitle::AbstractString = "", curves = [],
                    palette::Symbol = :abaqus, size = (1420, 820), frame::Int = nframes(v), field::Symbol = :vm_surf)
    return App(; title) do session
        WGLMakie.activate!()
        fd0 = v.fields[fieldindex(v, field)]
        fig = Figure(; size, backgroundcolor = :white, figure_padding = (4, 12, 8, 8))
        ls = LScene(fig[1, 1]; show_axis = false, scenekw = (; lights = _lights(), backgroundcolor = RGBf(0.96, 0.96, 0.95), clear = true))
        p = draw_model!(ls, v, frame, fd0; scale = v.default_scale, cmax = fd0.auto, palette)
        set_camera!(ls, v)
        gl = GridLayout(fig[1, 2])
        m1, m2 = charts!(gl, v, frame; curves)
        colsize!(fig.layout, 2, Fixed(420))

        nn = Base.size(v.x, 1); npv = length(v.pv_node)
        tohex(c) = "#" * hex(RGB{Colors.N0f8}(clamp(c.r, 0, 1), clamp(c.g, 0, 1), clamp(c.b, 0, 1)))
        data = Dict(
            "nn" => nn, "npv" => npv, "nf" => nframes(v), "nbands" => NBANDS,
            "x" => vec(permutedims(v.x)), "pv_node" => _i32(v.pv_node .- 1), "pv_off" => vec(permutedims(v.pv_off)),
            "sv_pv" => _i32(v.sv_pv .- 1), "sv_side" => Int8.(v.sv_side),
            "faces" => UInt32[UInt32(i - 1) for f in v.faces for i in f], "face_part" => _i32(v.face_part),
            "edges" => UInt32[UInt32(i - 1) for e in v.edges for i in e], "edge_part" => _i32(v.edge_part),
            "feature" => UInt32[UInt32(i - 1) for e in v.feature for i in e], "feature_part" => _i32(v.feature_part),
            "bolt_nodes" => _i32(v.bolt_nodes .- 1), "plate_t" => v.plate_t,
            "disp" => vec(v.disp),
            "stored" => Dict(String(k) => _u16(a) for (k, a) in v.stored),
            "range" => Dict(String(k) => [r[1], r[2]] for (k, r) in v.range),
            "fields" => [Dict("key" => String(f.key), "label" => f.label, "unit" => f.unit, "kind" => String(f.kind), "top" => String(f.top),
                              "bot" => String(f.bot), "auto" => f.auto, "plate_only" => f.plate_only) for f in v.fields],
            "seq" => tohex.(sequential_bands(palette)), "div" => tohex.(diverging_bands(palette)),
            "step" => v.step, "N" => v.N, "M" => v.M, "theta" => v.θ, "w" => v.w, "bolt" => vec(v.bolt), "nbolt" => Base.size(v.bolt, 1),
            "contact_max" => v.contact_max, "default_scale" => v.default_scale,
            "frame0" => frame - 1, "field0" => String(field),
        )

        sel = DOM.select([DOM.option(f.label; value = String(f.key)) for f in v.fields]...)
        slider = DOM.input(; type = "range", min = "0", max = string(nframes(v) - 1), step = "1", value = string(frame - 1))
        play = DOM.button("Play")
        steplab = DOM.span(""; style = "min-width: 70px; color: #0b0b0b;")
        scale = DOM.input(; type = "number", min = "0", step = "any", value = string(v.default_scale))
        cmax = DOM.input(; type = "number", min = "0", step = "any", value = string(fd0.auto))
        auto = DOM.button("Auto")
        cb(lbl, checked) = (i = DOM.input(; type = "checkbox", checked...); (i, DOM.label(i, lbl)))
        c_up, l_up = cb("Upright", (; checked = true)); c_ed, l_ed = cb("Mesh edges", (; checked = true))
        c_gh, l_gh = cb("Undeformed outline", (;)); c_bo, l_bo = cb("Bolts", (; checked = true))
        legend = DOM.div(; class = "legend"); info = DOM.div(; class = "info")

        code = js"""
        (() => {
        const D = $(data);
        const el = {sel: $(sel), slider: $(slider), play: $(play), steplab: $(steplab), scale: $(scale), cmax: $(cmax), auto: $(auto),
                    up: $(c_up), ed: $(c_ed), gh: $(c_gh), bo: $(c_bo), legend: $(legend), info: $(info)};
        const fieldByKey = {}; D.fields.forEach(f => fieldByKey[f.key] = f);
        el.sel.value = D.field0;
        // optional initial state from the URL hash: #frame=12&field=peeq&scale=10&cmax=60&upright=0&edges=0
        const hash = new URLSearchParams(window.location.hash.slice(1));
        if (hash.has('field') && fieldByKey[hash.get('field')]) el.sel.value = hash.get('field');
        if (hash.has('frame')) el.slider.value = Math.max(0, Math.min(D.nf - 1, parseInt(hash.get('frame'))));
        if (hash.has('scale')) el.scale.value = hash.get('scale');
        if (hash.has('upright')) el.up.checked = hash.get('upright') !== '0';
        if (hash.has('edges')) el.ed.checked = hash.get('edges') !== '0';
        const num = (x) => Number(Number(x).toPrecision(4));
        el.cmax.value = hash.has('cmax') ? hash.get('cmax') : num(fieldByKey[el.sel.value].auto);

        const nsv = D.sv_pv.length, NB = D.nbands;
        const pos = new Float32Array(nsv * 3), nor = new Float32Array(nsv * 3), col = new Float32Array(nsv);
        const pick = (idx, part, plateOnly, n) => {
            if (!plateOnly) return idx;
            const out = [];
            for (let i = 0; i < part.length; i++) if (part[i] === 1) for (let j = 0; j < n; j++) out.push(idx[n * i + j]);
            return new Uint32Array(out);
        };
        const fmt = (x) => { const a = Math.abs(x); if (x === 0) return '0';
            return a >= 100 ? x.toFixed(0) : a >= 10 ? x.toFixed(1) : a >= 1 ? x.toFixed(2) : Number(x.toPrecision(3)).toString(); };

        function positions(k, s) {
            const U = D.disp, o = k * D.nn * 3;
            for (let i = 0; i < nsv; i++) {
                const p = D.sv_pv[i], n = D.pv_node[p], sg = D.sv_side[i];
                pos[3*i]   = D.x[3*n]   + s * U[o + 3*n]   + sg * D.pv_off[3*p];
                pos[3*i+1] = D.x[3*n+1] + s * U[o + 3*n+1] + sg * D.pv_off[3*p+1];
                pos[3*i+2] = D.x[3*n+2] + s * U[o + 3*n+2] + sg * D.pv_off[3*p+2];
            }
        }
        function normals() {
            nor.fill(0);
            const F = D.faces;
            for (let f = 0; f < F.length; f += 3) {
                const a = 3*F[f], b = 3*F[f+1], c = 3*F[f+2];
                const ux = pos[b]-pos[a], uy = pos[b+1]-pos[a+1], uz = pos[b+2]-pos[a+2];
                const vx = pos[c]-pos[a], vy = pos[c+1]-pos[a+1], vz = pos[c+2]-pos[a+2];
                const nx = uy*vz-uz*vy, ny = uz*vx-ux*vz, nz = ux*vy-uy*vx;
                nor[a] += nx; nor[a+1] += ny; nor[a+2] += nz; nor[b] += nx; nor[b+1] += ny; nor[b+2] += nz; nor[c] += nx; nor[c+1] += ny; nor[c+2] += nz;
            }
            for (let i = 0; i < nsv; i++) {
                const l = Math.hypot(nor[3*i], nor[3*i+1], nor[3*i+2]) || 1;
                nor[3*i] /= l; nor[3*i+1] /= l; nor[3*i+2] /= l;
            }
        }
        function colors(k, f, cmax) {
            const lo = f.kind === 'seq' ? 0 : -cmax, hi = cmax, base = f.kind === 'seq' ? 0 : NB;
            const top = D.stored[f.top], bot = D.stored[f.bot], rt = D.range[f.top], rb = D.range[f.bot], o = k * D.npv;
            let vmax = -Infinity, vmin = Infinity;
            for (let i = 0; i < nsv; i++) {
                const t = D.sv_side[i] > 0, q = (t ? top : bot)[o + D.sv_pv[i]];
                if (q === 65535) { col[i] = NaN; continue; }
                const r = t ? rt : rb, val = r[0] + q / 65534 * (r[1] - r[0]);
                if (val > vmax) vmax = val; if (val < vmin) vmin = val;
                col[i] = base + NB * Math.max(0, Math.min(1 - 1e-6, (val - lo) / (hi - lo)));
            }
            return [vmin, vmax];
        }
        function segments(idx) {
            const out = new Float32Array(idx.length * 3);
            for (let i = 0; i < idx.length; i++) { const a = 3 * idx[i]; out[3*i] = pos[a]; out[3*i+1] = pos[a+1]; out[3*i+2] = pos[a+2]; }
            return out;
        }
        function legend(f, cmax, ext) {
            const cols = f.kind === 'seq' ? D.seq : D.div, lo = f.kind === 'seq' ? 0 : -cmax;
            let h = '<div class="ttl">' + f.label + (f.unit ? ' (' + f.unit + ')' : '') + '</div>';
            for (let b = NB - 1; b >= 0; b--)
                h += '<div class="row"><div class="sw" style="background:' + cols[b] + '"></div><div class="lab">' + fmt(lo + (b + 1) * (cmax - lo) / NB) + '</div></div>';
            h += '<div class="row last"><div style="width:30px"></div><div class="lab">' + fmt(lo) + '</div></div>';
            h += '<div class="note">this step: ' + (isFinite(ext[0]) ? fmt(ext[0]) + ' to ' + fmt(ext[1]) : 'no values') + '</div>';
            el.legend.innerHTML = h;
        }

        Promise.all([$(p.mesh), $(p.edges), $(p.feature), $(p.ghost), $(p.bolts), $(m1), $(m2)]).then(([pm, pe, pf, pg, pb, pk1, pk2]) => {
            const po = (a) => a[0].plot_object || a[0];   // find_plots returns the three.js objects
            const mesh = po(pm), edges = po(pe), feature = po(pf), ghost = po(pg), bolts = po(pb), mk1 = po(pk1), mk2 = po(pk2);
            const setbuf = (plot, names, arr) => {   // update the first attribute that exists
                const g = plot.mesh.geometry, key = names.find(n => (n in g.attributes) || (g.interleaved_attributes && n in g.interleaved_attributes));
                if (key) plot.update([[key, arr]]); else console.warn('BasePlateViewer: none of', names, 'in', Object.keys(g.attributes));
            };
            window.__bpv = {el, D, mesh, edges, feature, ghost, bolts, mk1, mk2};
            let lastUp = null;
            function render() {
                const k = parseInt(el.slider.value), f = fieldByKey[el.sel.value];
                const s = parseFloat(el.scale.value) || 0, cmax = parseFloat(el.cmax.value) > 0 ? parseFloat(el.cmax.value) : f.auto;
                const up = el.up.checked;
                positions(k, s); normals();
                const ext = colors(k, f, cmax);
                setbuf(mesh, ['positions_transformed_f32c', 'position'], pos);
                setbuf(mesh, ['normals', 'normal'], nor);
                setbuf(mesh, ['vertex_color', 'color'], col);
                if (lastUp !== up) {
                    mesh.update([['faces', pick(D.faces, D.face_part, !up, 3)]]);
                    lastUp = up;
                }
                edges.update([['positions_transformed_f32c', segments(pick(D.edges, D.edge_part, !up, 2))], ['visible', el.ed.checked]]);
                feature.update([['positions_transformed_f32c', segments(pick(D.feature, D.feature_part, !up, 2))]]);
                ghost.update([['visible', el.gh.checked]]);
                const bp = new Float32Array(D.bolt_nodes.length * 3), o = k * D.nn * 3;
                D.bolt_nodes.forEach((n, i) => { bp[3*i] = D.x[3*n] + s * D.disp[o+3*n]; bp[3*i+1] = D.x[3*n+1] + s * D.disp[o+3*n+1];
                                                 bp[3*i+2] = D.x[3*n+2] + s * D.disp[o+3*n+2] + D.plate_t / 2; });
                setbuf(bolts, ['positions_transformed_f32c', 'offset'], bp);
                bolts.update([['visible', el.bo.checked]]);
                setbuf(mk1, ['wgl_positions', 'positions_transformed_f32c'], new Float32Array([D.theta[k], D.M[k]]));
                setbuf(mk2, ['wgl_positions', 'positions_transformed_f32c'], new Float32Array([-D.w[k], D.N[k]]));
                legend(f, cmax, ext);
                el.steplab.textContent = 'Step ' + k + ' / ' + (D.nf - 1);
                let b = '';
                for (let i = 0; i < D.nbolt; i++) b += (i ? ', ' : '') + fmt(D.bolt[k * D.nbolt + i]);
                el.info.innerHTML = (D.step[k] === 1 ? 'Axial step' : 'Moment step') + ' &nbsp; <b>N</b> ' + fmt(D.N[k]) + ' kip &nbsp; <b>M</b> ' + fmt(D.M[k]) +
                    ' kip-in &nbsp; <b>θ</b> ' + D.theta[k].toFixed(5) + ' rad &nbsp; <b>bolt forces</b> ' + b + ' kip &nbsp; <b>max contact</b> ' + fmt(D.contact_max[k]) + ' ksi';
                window.__bpv_rendered = (window.__bpv_rendered || 0) + 1;
            }
            let timer = null;
            const stop = () => { if (timer) { clearInterval(timer); timer = null; el.play.textContent = 'Play'; } };
            el.play.onclick = () => {
                if (timer) return stop();
                if (parseInt(el.slider.value) >= D.nf - 1) el.slider.value = 0;
                el.play.textContent = 'Pause';
                timer = setInterval(() => { const k = parseInt(el.slider.value); if (k >= D.nf - 1) return stop(); el.slider.value = k + 1; render(); }, 180);
            };
            el.slider.oninput = render; el.scale.oninput = render; el.cmax.oninput = render;
            el.sel.onchange = () => { el.cmax.value = num(fieldByKey[el.sel.value].auto); render(); };
            el.auto.onclick = () => { el.cmax.value = num(fieldByKey[el.sel.value].auto); render(); };
            [el.up, el.ed, el.gh, el.bo].forEach(c => c.onchange = render);
            render();
        }).catch(e => console.error('BasePlateViewer:', e));
        })();
        """
        Bonito.onload(session, legend, js"() => { $(code) }")
        return DOM.div(DOM.style(CSS),
            DOM.h1(title), DOM.div(subtitle; class = "sub"),
            DOM.div(DOM.label("Contour ", sel), DOM.label("max ", cmax, auto), DOM.label(play, slider, steplab; class = "step"),
                    DOM.label("Deformation × ", scale), l_up, l_ed, l_gh, l_bo; class = "bar"),
            DOM.div(DOM.div(legend; class = "left"), DOM.div(fig, info; class = "right"); class = "wrap"); class = "bpv")
    end
end
viewer_app(res::BasePlateFEA.Result; max_frames::Int = 40, second_order::Bool = false, kw...) =
    viewer_app(VizData(res; max_frames, second_order); kw...)

"""
    serve_viewer(result; port = 8642, kwargs...) -> Bonito.Server

Serve the viewer at `http://localhost:port` (keywords as [`viewer_app`](@ref)).
"""
function serve_viewer(res; port::Int = 8642, kw...)
    app = viewer_app(res; kw...)
    server = Bonito.Server(app, "127.0.0.1", port)
    println("BasePlateViewer: http://localhost:$port")
    return server
end

"""
    export_viewer(path, result; kwargs...) -> path

Write the viewer as one standalone HTML file (keywords as [`viewer_app`](@ref)). The page needs no
Julia process and no internet connection; all interaction runs in the browser.
"""
function export_viewer(path::AbstractString, res; kw...)
    Bonito.export_static(path, viewer_app(res; kw...))
    return path
end

end # module

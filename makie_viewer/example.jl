# Example: viewer page and still images from a saved result.
#   julia --project=makie_viewer makie_viewer/example.jl <result.jls> <output dir> [serve]
# <result.jls> is a `BasePlateFEA.Result` written with `Serialization.serialize`.
using Serialization, BasePlateFEA
include(joinpath(@__DIR__, "BasePlateViewer.jl"))
using .BasePlateViewer

res = deserialize(ARGS[1])
out = mkpath(ARGS[2])
name = splitext(basename(ARGS[1]))[1]
v = VizData(res; max_frames = 40, second_order = false)          # precompute once, reuse below

html = export_viewer(joinpath(out, name * ".html"), v; title = name)
println(html, "  ", round(filesize(html) / 1e6, digits = 1), " MB")
println(render_image(joinpath(out, name * "_vm.png"), v; field = :vm_surf, title = name, charts = true))
println(render_image(joinpath(out, name * "_peeq.png"), v; field = :peeq, title = name))

if length(ARGS) > 2 && ARGS[3] == "serve"
    server = serve_viewer(v; title = name)
    wait(server)
end

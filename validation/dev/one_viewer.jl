# one viewer page:  julia --project=makie_viewer validation/dev/one_viewer.jl <vdir> <resdir> <name> <second_order> <out.html>
using Serialization, DelimitedFiles, BasePlateFEA
include(joinpath(@__DIR__, "..", "..", "makie_viewer", "BasePlateViewer.jl")); using .BasePlateViewer
res = deserialize(joinpath(ARGS[2], ARGS[3] * ".jls"))
v = VizData(res; max_frames = 40, second_order = ARGS[4] == "true")
export_viewer(ARGS[5], v; title = ARGS[3]); println("ok ", filesize(ARGS[5]))

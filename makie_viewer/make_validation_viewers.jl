# Makie viewer pages for the Barcelona validation runs, with the test curves overlaid.
#   julia --project=makie_viewer makie_viewer/make_validation_viewers.jl <validation dir> <results dir> <suffix> [second_order] [max_frames]
# Every `<type>_N<load>_<suffix>.jls` in <results dir> becomes `<validation dir>/viewers_makie/<same name>.html`.
using Serialization, DelimitedFiles, BasePlateFEA
include(joinpath(@__DIR__, "BasePlateViewer.jl"))
using .BasePlateViewer

function test_curves(vdir, typ, N)
    datadir = joinpath(vdir, "barcelona_test_data"); cs = []
    for f in sort(readdir(datadir))
        startswith(f, "type$(typ)_P$(round(Int, 1000N))_") && endswith(f, ".csv") || continue
        d = readdlm(joinpath(datadir, f), ',', Float64; skipstart = 1); n = argmax(d[:, 1])
        push!(cs, (label = "test", kind = :test, θ = d[1:n, 1], M = d[1:n, 2]))
    end
    return cs
end

"every computed step from the run's CSV (the saved result may hold fewer frames)"
function full_history(csv)
    isfile(csv) || return nothing
    hdr = split(readline(csv), ','); d = readdlm(csv, ',', Float64; skipstart = 1)
    col(n) = d[:, findfirst(==(n), hdr)]
    return (step = round.(Int, col("step")), N = col("N_kip"), M = col("M_kipin"), θ = col("rotation_rad"), w = col("w_in"))
end

function main(vdir, resdir, suffix, second_order = false, max_frames = 40)
    out = mkpath(joinpath(vdir, "viewers_makie"))
    for f in sort(readdir(resdir))
        m = match(Regex("^([ABC])_N([0-9p]+)_$(suffix)\\.jls\$"), f)
        m === nothing && continue
        typ = m[1]; N = parse(Float64, replace(m[2], "p" => "."))
        name = splitext(f)[1]; path = joinpath(out, name * ".html")
        res = deserialize(joinpath(resdir, f))
        v = VizData(res; max_frames, second_order, history = full_history(joinpath(resdir, name * ".csv")),
                    actuator_length = 1250 / 25.4)   # axial jack of the Barcelona rig, hinge to hinge, 1250 mm
        export_viewer(path, v; title = "Type $typ, N = $(m[2] |> x -> replace(x, "p" => ".")) kip",
                      subtitle = name, curves = test_curves(vdir, typ, N))
        println(name, ".html  ", round(filesize(path) / 1e6, digits = 1), " MB"); flush(stdout)
        res = v = nothing; GC.gc()
    end
end

main(ARGS[1], ARGS[2], ARGS[3], length(ARGS) > 3 && ARGS[4] == "true", length(ARGS) > 4 ? parse(Int, ARGS[5]) : 40)

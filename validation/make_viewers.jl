# Write the 3D viewer page of every saved Barcelona result, with the test curves overlaid.
#   julia --project=. validation/make_viewers.jl <validation dir> [name filter]
using BasePlateFEA, Serialization, DelimitedFiles
include(joinpath(@__DIR__, "barcelona_cases.jl"))

"test curves of a plate type and nominal axial load (kip), trimmed at the largest rotation"
function test_curves(datadir, typ, N)
    out = []
    for f in sort(readdir(datadir))
        startswith(f, "type$(typ)_P$(round(Int, 1000N))_") && endswith(f, ".csv") || continue
        d = readdlm(joinpath(datadir, f), ',', Float64; skipstart = 1)
        n = argmax(d[:, 1])
        push!(out, (label = replace(f, ".csv" => ""), kind = :test, θ = d[1:n, 1], M = d[1:n, 2]))
    end
    return out
end

function main(vdir, filter = "")
    resdir = joinpath(vdir, "results"); outdir = mkpath(joinpath(vdir, "viewers"))
    for f in sort(readdir(resdir))
        endswith(f, ".jls") && occursin(filter, f) || continue
        name = replace(f, ".jls" => ""); typ, nstr, mode = split(name, "_")
        N = parse(Float64, replace(nstr[2:end], "p" => "."))
        res = deserialize(joinpath(resdir, f))
        curves = test_curves(joinpath(vdir, "barcelona_test_data"), typ, N)
        if mode == "plastic" && isfile(joinpath(resdir, "$(typ)_$(nstr)_elastic.csv"))
            e = readdlm(joinpath(resdir, "$(typ)_$(nstr)_elastic.csv"), ',', Float64; skipstart = 1)
            e = e[findlast(==(1), e[:, 1]):end, :]
            push!(curves, (label = "Elastic model", kind = :reference, θ = e[:, 4], M = e[:, 3]))
        end
        c = BARCELONA[Symbol(typ)]
        sub = mode == "plastic" ? "Axial load $N kip · plate Fy = $(round(c.plate[1], digits = 1)) ksi, upright Fy = $(round(c.upright[1], digits = 1)) ksi · $(res.status)" :
                                  "Axial load $N kip · elastic steel (bp_Barcelona.m settings) · $(res.status)"
        p = write_viewer(joinpath(outdir, name * ".html"), res; title = "Barcelona base plate type $typ", subtitle = sub, curves)
        println(p, "  ", round(filesize(p) / 1e6, digits = 1), " MB")
    end
end

main(ARGS...)

# Initial-stiffness sensitivity: one Barcelona case to θ = 0.012 with one support parameter changed.
#   julia --project=. validation/run_stiffness_sensitivity.jl <outdir> <A|B|C> <N kip> <tag> [key=value ...]
# keys: Ec (ksi), contact_stiffness (ksi/in), bolt_stiffness (kip/in)
using BasePlateFEA, LinearAlgebra, Printf
include(joinpath(@__DIR__, "barcelona_cases.jl"))
BLAS.set_num_threads(1)
outdir, typ, N, tag = ARGS[1], Symbol(ARGS[2]), parse(Float64, ARGS[3]), ARGS[4]
kw = Dict(Symbol(split(a, "=")[1]) => parse(Float64, split(a, "=")[2]) for a in ARGS[5:end])
res = analyze(barcelona_input(typ; plastic = true, kw...); N, Mmax = 500.0, θmax = 0.012, dM0 = 0.2, dM_max = 2.0)
open(joinpath(outdir, "$(typ)_N$(loadtag(N))_$(tag).csv"), "w") do io
    println(io, "step,M_kipin,rotation_rad,bolt1_kip,bolt2_kip")
    for f in res.frames
        @printf(io, "%d,%.6g,%.8g,%.6g,%.6g\n", f.step, f.M, f.θ, f.bolt_force[1], f.bolt_force[2])
    end
end
println("$(typ) $(tag): $(res.status), $(length(res.frames)) frames")

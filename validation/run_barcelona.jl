# Run one Barcelona case and save the curve (CSV) and the full result (serialized).
#   julia --project=. validation/run_barcelona.jl <outdir> <A|B|C> <N kip> <elastic|plastic> [pure]
# With `pure` the moment is applied without the lateral force (Ht → ∞); the file name gets `_pure`.
using BasePlateFEA, Serialization, LinearAlgebra, Printf
include(joinpath(@__DIR__, "barcelona_cases.jl"))
BLAS.set_num_threads(2)

function save_curve(path, res)
    open(path, "w") do io
        println(io, "step,N_kip,M_kipin,rotation_rad,w_in,bolt1_kip,bolt2_kip,max_contact_ksi,max_peeq,max_vm_plate_ksi,iterations")
        m = res.model.mesh
        plate = findall(==(1), m.part)
        for f in res.frames
            mx(v) = maximum(x -> isnan(x) ? 0.0 : x, v)
            vm = max(mx(f.fields[:vm_top][plate]), mx(f.fields[:vm_bot][plate]))
            @printf(io, "%d,%.6g,%.6g,%.8g,%.8g,%.6g,%.6g,%.6g,%.6g,%.6g,%d\n", f.step, f.N, f.M, f.θ, f.w,
                    f.bolt_force[1], f.bolt_force[2], maximum(f.contact_stress), mx(f.fields[:peeq]), vm, f.iterations)
        end
    end
end

function main(outdir, typ, N, mode, pure = false)
    plastic = mode == "plastic"
    inp = pure ? barcelona_input(typ; plastic, Ht = 1e9) : barcelona_input(typ; plastic)
    t0 = time()
    res = plastic ? analyze(inp; N, Mmax = 500.0, θmax = BARCELONA[typ].θmax, dM0 = 0.2, dM_max = 5.0) :
                    analyze(inp; N, Mmax = barcelona_Mmax0(typ) + N)
    name = "$(typ)_N$(loadtag(N))_$(mode)" * (pure ? "_pure" : "")
    save_curve(joinpath(outdir, name * ".csv"), res)
    serialize(joinpath(outdir, name * ".jls"), res)
    @printf("%s: %s, %d frames, %.0f s\n", name, res.status, length(res.frames), time() - t0)
end

main(ARGS[1], Symbol(ARGS[2]), parse(Float64, ARGS[3]), ARGS[4], length(ARGS) > 4 && ARGS[5] == "pure")

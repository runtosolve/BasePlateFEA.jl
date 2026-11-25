# Continue a stopped tall run with another control:  julia --project=. validation/dev/resume_arclength.jl <in.jls> <out.jls> <θmax> <control> [steps]
using Serialization, BasePlateFEA
res = deserialize(ARGS[1]); θmax = parse(Float64, ARGS[3])
keep = findlast(f -> f.θ < 1 && abs(f.M) < 1e4, res.frames)        # drop a runaway last frame (its state is not the committed one, so only the log is meaningful then)
keep == length(res.frames) || println("warning: last frame is a runaway state; shell state does not match frame $keep")
r0 = BasePlateFEA.Result(res.model, res.frames[1:keep], res.status)
out = analyze(res.model; N = r0.frames[end].N, Mmax = 500.0, θmax, dM0 = 0.2, dM_max = 2.0, dθ_max = θmax / 80,
              max_steps = length(ARGS) > 4 ? parse(Int, ARGS[5]) : 150, control = Symbol(ARGS[4]), resume = r0, linesearch = get(ENV, "BPFEA_LS", "") == "1", trace_failures = get(ENV, "BPFEA_TRACE", "") == "1")
println(out.status, "  frames ", length(out.frames), "  θ = ", out.frames[end].θ)
serialize(ARGS[2], out)

# Type A with the upright modelled to the load point (items 1–3 of the Nov 22 review):
#   julia --project=. validation/run_typeA_tall.jl <outdir> <N kip> <tag> [imperfection/t] [bolt stiffness kip/in] [θmax] [type A|B|C] [max steps] [pure]
# With `pure` the moment is applied at the loaded end as a couple, without the lateral force (Ht → ∞).
# Upright 24 in tall (fine mesh over the lower 9 in, elastic above), loaded end held in the plane
# of bending, rotation measured 85 mm above the plate, Yun–Gardner steel curves, imperfection in
# the shape of the first buckling mode, base moment including N·δ.
using BasePlateFEA, Serialization, LinearAlgebra, Printf
include(joinpath(@__DIR__, "barcelona_cases.jl"))
BLAS.set_num_threads(1)

function tall_input(typ; bolt_stiffness = 120.0, kwargs...)
    c = BARCELONA[typ]
    return BasePlateInput(; geometry = barcelona_test_geometry(typ),
                          plate_steel = yun_gardner_steel(c.plate[1], c.plate[2]),
                          upright_steel = yun_gardner_steel(c.upright[1], c.upright[2]),
                          H = 24.0, H_fine = 9.0, fe_size_upper = 0.5, upright_plastic_height = 9.0,
                          top_restraint = :planar, Hm = 85 / 25.4, bolt_stiffness, kwargs...)
end

function save_curve(path, res)
    open(path, "w") do io
        println(io, "step,N_kip,M_first_order_kipin,M_kipin,rotation_rad,w_in,top_displacement_in,bolt1_kip,bolt2_kip,max_contact_ksi,max_peeq_upright,max_peeq_plate,iterations")
        m = res.model.mesh
        plate = findall(==(1), m.part); upr = findall(>(1), m.part)
        mx(v) = maximum(x -> isnan(x) ? 0.0 : x, v)
        for f in res.frames
            @printf(io, "%d,%.6g,%.6g,%.6g,%.8g,%.8g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%d\n", f.step, f.N, f.M, base_moment(res.model, f),
                    f.θ, f.w, f.u[6*(m.master-1)+2], f.bolt_force[1], f.bolt_force[2], maximum(f.contact_stress),
                    mx(f.fields[:peeq][upr]), mx(f.fields[:peeq][plate]), f.iterations)
        end
    end
end

"result with at most `n` frames, spaced evenly in rotation (the unloaded state, the end of the axial step and the last step are kept)"
function thin(res, n = 45)
    nf = length(res.frames); i1 = findlast(f -> f.step == 1, res.frames)
    θ = [f.θ for f in res.frames[i1:end]]; s = cumsum(vcat(0.0, abs.(diff(θ))))
    pick = [i1 - 1 + searchsortedfirst(s, x) for x in range(0, s[end]; length = n)]
    keep = sort!(unique(vcat(1, i1, nf, clamp.(pick, i1, nf))))
    return Result(res.model, res.frames[keep], res.status)
end

function main(outdir, N, tag, imp = 0.0, kb = 120.0, θmax = 0.08, typ = :A, max_steps = 400, pure = false)
    # optional overrides: BPFEA_BOLT_FY (bolt yield, kip), BPFEA_BOLT_HARD (ultimate/yield), BPFEA_MODE (buckling mode used as imperfection)
    extra = (; (haskey(ENV, "BPFEA_BOLT_FY") ? (; bolt_yield = parse(Float64, ENV["BPFEA_BOLT_FY"]), bolt_hardening = parse(Float64, get(ENV, "BPFEA_BOLT_HARD", "1.25"))) : (;))...)
    if haskey(ENV, "BPFEA_BOLT_LAW")   # anchor law as "d1,f1;d2,f2;..." (in, kip) for the points after (0, 0)
        pts = [parse.(Float64, split(p, ',')) for p in split(ENV["BPFEA_BOLT_LAW"], ';')]
        extra = (; extra..., bolt_law = SpringLaw(vcat(-1000.0, 0.0, first.(pts)), vcat(0.0, 0.0, last.(pts))))
    end
    mode = parse(Int, get(ENV, "BPFEA_MODE", "1"))
    inp = pure ? tall_input(typ; bolt_stiffness = kb, Ht = 1e9, extra...) : tall_input(typ; bolt_stiffness = kb, extra...)
    mesh = generate_mesh(inp)
    t0 = time()
    if imp != 0
        λ, Φ = buckling_modes(Model(inp; mesh); N, M = 40.0, nmodes = max(4, mode))
        @printf("buckling load factors for N = %g kip, M = 40 kip-in: %s\n", N, string(round.(λ, digits = 3)))
        mesh = apply_imperfection(mesh, Φ[:, mode], imp * inp.geometry.tu)
    end
    res = analyze(Model(inp; mesh); N, Mmax = 500.0, θmax, dM0 = 0.2, dM_max = 2.0, dθ_max = θmax / 80, max_steps, trace_failures = get(ENV, "BPFEA_TRACE", "") == "1", control = Symbol(get(ENV, "BPFEA_CONTROL", "displacement")))
    name = "$(typ)_N$(loadtag(N))_$(tag)"
    save_curve(joinpath(outdir, name * ".csv"), res)
    serialize(joinpath(outdir, name * ".jls"), thin(res))
    @printf("%s: %s, %d frames, %.0f s\n", name, res.status, length(res.frames), time() - t0)
end

main(ARGS[1], parse(Float64, ARGS[2]), ARGS[3], (parse(Float64, a) for a in ARGS[4:min(end, 6)])...,
     (length(ARGS) > 6 ? (Symbol(ARGS[7]),) : ())..., (length(ARGS) > 7 ? (parse(Int, ARGS[8]),) : ())..., (length(ARGS) > 8 ? (ARGS[9] == "pure",) : ())...)

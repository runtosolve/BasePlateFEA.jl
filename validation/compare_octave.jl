# Compare BasePlateFEA with S. Ádány's MATLAB code (run in GNU Octave) on the identical mesh:
# stiffness matrices, axial shortening and the moment–rotation curve of every reference case.
#   julia --project=. validation/compare_octave.jl <validation dir>
using BasePlateFEA, SparseArrays, LinearAlgebra, DelimitedFiles, Printf
const B = BasePlateFEA
interp(x, xs, ys) = (i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1); ys[i] + (ys[i+1] - ys[i]) * (x - xs[i]) / (xs[i+1] - xs[i]))

function main(vdir)
    dir = joinpath(vdir, "octave_reference"); out = mkpath(joinpath(vdir, "results"))
    io = open(joinpath(out, "octave_comparison.csv"), "w")
    println(io, "case,dK_steel,dK_beam,dK_support,dK_concrete,w_axial_ratio,K0_octave,K0_julia,theta_ratio_at_half_Mmax,theta_octave_at_Mmax,theta_julia_at_Mmax,theta_octave_after_one_newton_correction,octave_residual_over_step_load")
    for typ in (:A, :B, :C), N in (5.0, 10.0, 15.0, 20.0)
        tag = "$(typ)_N$(Int(N))"
        inp = BasePlateInput(geometry = barcelona_geometry(typ))
        mesh = mesh_from_matlab(dir, tag, inp)
        nn = B.nnodes(mesh); ncb = length(mesh.bp_nodes)
        rdK(name, n) = (t = readdlm(joinpath(dir, "$(tag)_$(name).txt"), ','); sparse(Int.(t[:, 1]), Int.(t[:, 2]), t[:, 3], n, n))
        rel(A, R) = norm(A - R) / norm(R)
        model = Model(inp; mesh)
        dK = (rel(model.shell.K_lin, rdK("Ke_steel", 6nn)), rel(B.beam_stiffness(mesh, inp), rdK("Ke_beam", 6nn)),
              rel(B.support_stiffness(mesh, inp), rdK("Ke_sup", 6nn)), rel(B.concrete_stiffness(mesh, inp), rdK("Ke_concr", ncb)))
        # equilibrium of the Octave end state, evaluated in the Julia model
        Mmax = barcelona_Mmax0(typ) + N
        H, Ht, _ = B.rounded_heights(inp)
        P1 = B.master_load(mesh, (0, 0, -N, 0, 0, 0)); P2 = B.master_load(mesh, (0, 1 / Ht, 0, -(1 - H / Ht), 0, 0))
        d = vec(readdlm(joinpath(dir, "$(tag)_d_final.txt"), ','))
        uo = vcat(d[ncb+1:end], d[1:ncb])
        f, K = B.internal_forces(model, uo)
        r = f - P1 - Mmax * P2
        θcorr = abs(B.rotation(model, uo - B.factorize_tangent(K) \ r)[1])
        res = analyze(model; N, Mmax, verbose = false)
        θ, M = moment_rotation(res)
        ref = readdlm(joinpath(dir, "$(tag)_m_r.txt"), ','); ref = ref[ref[:, 1] .> 0, :]
        nd = readdlm(joinpath(dir, "$(tag)_n_d_step1.txt"), ',')
        w1 = res.frames[findlast(f -> f.step == 1, res.frames)].w
        writedlm(joinpath(out, "$(tag)_julia_on_octave_mesh.csv"), vcat(["rotation_rad" "M_kipin"], hcat(θ, M)), ',')
        k = findfirst(>=(Mmax / 2), ref[:, 1])
        @printf(io, "%s,%.2e,%.2e,%.2e,%.2e,%.5f,%.1f,%.1f,%.4f,%.6f,%.6f,%.6f,%.4f\n", tag, dK..., w1 / nd[end, 2],
                ref[2, 1] / ref[2, 2], interp(ref[2, 1], M, θ) |> t -> ref[2, 1] / t, interp(ref[k, 1], M, θ) / ref[k, 2],
                ref[end, 2], θ[end], θcorr, norm(r) / norm(Mmax * P2))
        flush(io)
        println(tag, " done: ", res.status)
    end
    close(io)
end
main(ARGS[1])

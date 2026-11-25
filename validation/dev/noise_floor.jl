# Is the Newton stall a round-off floor?  julia --project=. validation/dev/noise_floor.jl <result.jls>
using Serialization, LinearAlgebra, SparseArrays, Printf, BasePlateFEA
const B = BasePlateFEA
res = deserialize(ARGS[1]); model = res.model; m = model.mesh; fr = res.frames[end]
u = copy(fr.u); n = length(u)
H, Ht, _ = B.rounded_heights(model.inp)
P1 = B.master_load(m, (0.0, 0.0, -fr.N, 0.0, 0.0, 0.0)); P2 = B.master_load(m, (0.0, 1 / Ht, 0.0, -(1 - H / Ht), 0.0, 0.0))
resid(u, λ) = (B.restore!(model.shell); (f, K) = B.internal_forces(model, u); (f .- P1 .- λ .* P2, f, K))
r, f, K = resid(u, fr.M); fref = max(norm(P1 .+ fr.M .* P2), norm(f))
@printf("θ = %.5f  M = %.4f  |r|/fref = %.3e  fref = %.4g  |u|max = %.3g\n", fr.θ, fr.M, norm(r) / fref, fref, maximum(abs, u))
r2, _, _ = resid(u, fr.M); @printf("repeatability |r2-r|/fref = %.3e\n", norm(r2 - r) / fref)
Kc = model.K_const
@printf("round-off estimate eps*| |Kc| |u| | / fref = %.3e   max diag Kc = %.3e  max diag K = %.3e\n", eps() * norm(abs.(Kc) * abs.(u)) / fref, maximum(diag(Kc)), maximum(diag(K)))
p = sortperm(abs.(r), rev = true)[1:8]
for i in p
    i <= 6 * B.nnodes(m) ? @printf("  dof %d node %d dir %d r = %.3e  x = %s%s\n", i, (i - 1) ÷ 6 + 1, (i - 1) % 6 + 1, r[i], string(round.(size(m.x, 2) == 3 ? m.x[(i-1)÷6+1, :] : m.x[:, (i-1)÷6+1], digits = 3)), (i - 1) ÷ 6 + 1 == m.master ? " (master)" : "") :
                           @printf("  concrete dof %d r = %.3e\n", i, r[i])
end
# Newton iterations at fixed load: does the residual go below the floor?
c = P2; q = dot(c, u); λ = fr.M
for it in 1:8
    global u, λ
    r, f, K = resid(u, λ); F = B.factorize_tangent(K)
    dr = F \ r; dl = F \ P2; dλ = (q - dot(c, u) + dot(c, dr)) / dot(c, dl)
    du = dλ .* dl .- dr
    @printf("it %d |r|/fref = %.3e  |du|/|u| = %.3e  dλ = %.3e  linear-solve residual |K du + r - dλ P2|/fref = %.3e\n", it, norm(r) / fref, norm(du) / norm(u), dλ, norm(K * du .+ r .- dλ .* P2) / fref)
    u .+= du; λ += dλ
end

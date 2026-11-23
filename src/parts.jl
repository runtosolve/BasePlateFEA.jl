# --------------------------------------------------------------------------------------------
# Linear parts of the model: rigid-link beams (ke_beam_EB.m, stiffmat_e_beam.m), in-plane support
# springs (fem_support.m), the concrete half-space (stiffmat_concr.m), and the one-sided nonlinear
# springs of the anchor bolts and the plate–concrete contact links (spr_stif.m, spr_force.m,
# spr_update.m, spring_support.m, stifmat_cond.m, force_back.m).
#
# Global dofs: 6 per node, (node − 1)·6 + [ux uy uz θx θy θz], followed by one vertical concrete
# dof per base plate node (in the order of `mesh.bp_nodes`).
# --------------------------------------------------------------------------------------------

ndofs_steel(m::Mesh) = 6 * nnodes(m)
ndofs(m::Mesh) = 6 * nnodes(m) + length(m.bp_nodes)
nodedofs(i::Int) = (6 * (i - 1) + 1):(6 * i)

"""
    ke_beam_EB(P1, P2, A, Iy, Iz, J, E, G)

12×12 global stiffness matrix of the Euler–Bernoulli beam used for the rigid links, as in
`ke_beam_EB.m` (local x along the beam, local z as close as possible to global Z). Note that the
MATLAB code uses the torsional stiffness `G J / L³` (not `G J / L`); this is reproduced.
"""
function ke_beam_EB(P1, P2, A, Iy, Iz, J, E, G)
    L = norm(P2 - P1)
    j1 = (P2 - P1) / L
    j3 = cross(P2 - P1, [0.0, 0.0, 100.0])
    j3 /= norm(j3)
    j2 = cross(j3, j1)
    T = hcat(j1, j2, j3)
    kx = E * A / L; kxx = G * J / L^3
    ky2 = 2E * Iy / L; ky4 = 4E * Iy / L; ky6 = 6E * Iy / L^2; ky12 = 12E * Iy / L^3
    kz2 = 2E * Iz / L; kz4 = 4E * Iz / L; kz6 = 6E * Iz / L^2; kz12 = 12E * Iz / L^3
    ke = [
        kx 0 0 0 0 0 -kx 0 0 0 0 0
        0 kz12 0 0 0 kz6 0 -kz12 0 0 0 kz6
        0 0 ky12 0 -ky6 0 0 0 -ky12 0 -ky6 0
        0 0 0 kxx 0 0 0 0 0 -kxx 0 0
        0 0 -ky6 0 ky4 0 0 0 ky6 0 ky2 0
        0 kz6 0 0 0 kz4 0 -kz6 0 0 0 kz2
        -kx 0 0 0 0 0 kx 0 0 0 0 0
        0 -kz12 0 0 0 -kz6 0 kz12 0 0 0 -kz6
        0 0 -ky12 0 ky6 0 0 0 ky12 0 ky6 0
        0 0 0 -kxx 0 0 0 0 0 kxx 0 0
        0 0 -ky6 0 ky2 0 0 0 ky6 0 ky4 0
        0 kz6 0 0 0 kz2 0 -kz6 0 0 0 kz4
    ]
    TT = zeros(12, 12)
    for k in 0:3
        TT[3k+1:3k+3, 3k+1:3k+3] = T
    end
    return TT * ke * TT'
end

"global stiffness matrix of the rigid-link beams (steel dofs)"
function beam_stiffness(m::Mesh, inp::BasePlateInput)
    E = inp.upright_steel.E; G = E / 2 / (1 + inp.upright_steel.ν)
    I = Int[]; J = Int[]; V = Float64[]
    for (n1, n2, kind) in m.beams
        s = kind == 1 ? inp.top_mpc : inp.hole_mpc
        ke = ke_beam_EB(m.x[n1, :], m.x[n2, :], s.A, s.Iy, s.Iz, s.J, E, G)
        dofs = vcat(nodedofs(n1), nodedofs(n2))
        for b in 1:12, a in 1:12
            ke[a, b] == 0 && continue
            push!(I, dofs[a]); push!(J, dofs[b]); push!(V, ke[a, b])
        end
    end
    n = ndofs_steel(m)
    return sparse(I, J, V, n, n)
end

"""
    support_stiffness(m, inp)

Diagonal in-plane springs (global x and y) of `SUPPORT_table`: `hole_stiffness` on the hole
perimeter nodes and the 'friction' `friction_stiffness` on every base plate node.
"""
function support_stiffness(m::Mesh, inp::BasePlateInput)
    k = zeros(ndofs_steel(m))
    for i in m.bp_nodes
        k[6(i-1)+1] += inp.friction_stiffness
        k[6(i-1)+2] += inp.friction_stiffness
    end
    for hn in m.hole_nodes, i in hn
        k[6(i-1)+1] += inp.hole_stiffness
        k[6(i-1)+2] += inp.hole_stiffness
    end
    if inp.top_restraint == :planar                 # master node: ux, θy, θz held
        for d in (1, 5, 6)
            k[6(m.master-1)+d] += 1e9
        end
    elseif inp.top_restraint != :free
        throw(ArgumentError("top_restraint must be :free or :planar"))
    end
    return spdiagm(0 => k)
end

"""
    concrete_stiffness(m, inp) -> SparseMatrixCSC

Stiffness matrix of the concrete block on its vertical dofs: the inverse of the Boussinesq
flexibility matrix of an elastic half-space, where the unit force of base plate node `i` is spread
over a square patch of side `patch_factor × √(area_i)` centred on the node (settlement under and
around a uniformly loaded rectangle, Love 1929). Entries smaller than `concrete_threshold` times
the largest entry are dropped, as in `stiffmat_concr.m`.
"""
function concrete_stiffness(m::Mesh, inp::BasePlateInput)
    n = length(m.bp_nodes)
    vx = m.x[m.bp_nodes, 1]; vy = m.x[m.bp_nodes, 2]
    F = zeros(n, n)
    coef = (1 - inp.νc^2) / inp.Ec / π
    for ro in 1:n
        a = sqrt(m.bp_area[ro]) * inp.patch_factor / 2; b = a
        p = 1 / (4 * a * b)
        for co in 1:n
            x = vx[co] - vx[ro]; y = vy[co] - vy[ro]
            uu = (x + a) * log((y + b + hypot(y + b, x + a)) / (y - b + hypot(y - b, x + a)))
            uu += (y + b) * log((x + a + hypot(y + b, x + a)) / (x - a + hypot(y + b, x - a)))
            uu += (x - a) * log((y - b + hypot(y - b, x - a)) / (y + b + hypot(y + b, x - a)))
            uu += (y - b) * log((x - a + hypot(y - b, x - a)) / (x + a + hypot(y - b, x + a)))
            F[ro, co] = uu * p * coef
        end
    end
    K = inv(F)
    if inp.concrete_threshold > 0
        K[abs.(K) .< maximum(K) * inp.concrete_threshold] .= 0.0
    end
    return sparse(K)
end

# ---- one-sided springs with unloading ------------------------------------------------------

"force and tangent stiffness of the virgin curve at displacement `d`"
function virgin(law::SpringLaw, d)
    d <= law.d[2] && return 0.0, 0.0
    n = length(law.d)
    d >= law.d[n] && return law.f[n], 0.0
    i = 2
    while d >= law.d[i+1]
        i += 1
    end
    k = (law.f[i+1] - law.f[i]) / (law.d[i+1] - law.d[i])
    return law.f[i] + k * (d - law.d[i]), k
end

"""
    spring_response(law, d, dmax) -> (f, k)

Force and tangent stiffness at displacement `d` of a spring whose largest displacement so far is
`dmax`: the virgin curve for `d ≥ dmax`, otherwise elastic unloading with the initial stiffness
from the point reached at `dmax`, and no force below the displacement where it reaches zero.
"""
function spring_response(law::SpringLaw, d, dmax)
    (dmax <= 0 || d >= dmax) && return virgin(law, d)
    fmax, _ = virgin(law, dmax)
    d0 = dmax - fmax / law.k0
    return d <= d0 ? (0.0, 0.0) : (law.k0 * (d - d0), law.k0)
end

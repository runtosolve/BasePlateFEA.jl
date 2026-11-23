# --------------------------------------------------------------------------------------------
# Model assembly and the two-step incremental analysis (bp_Barcelona.m, FEanal_incr.m):
# step 1 applies the axial force, step 2 the lateral force and moment at the master node.
#
# Unlike the MATLAB code (one damped modified-Newton pass per increment, load control), every
# increment is iterated to equilibrium with the consistent tangent (shell plasticity, spring and
# contact states), and step 2 is driven by the work-conjugate displacement of the load pattern,
# so the analysis can follow the moment–rotation curve past the onset of a plateau.
# --------------------------------------------------------------------------------------------

"""
    Model(inp::BasePlateInput; mesh = generate_mesh(inp))

Assembled model: mesh, shell elements, the constant part of the stiffness matrix (linear shell
elements, rigid links, in-plane supports, concrete), the spring laws, and the spring history
(largest displacement of every bolt spring and contact link).
"""
struct Model
    inp::BasePlateInput
    mesh::Mesh
    shell::ShellModel
    K_const::SparseMatrixCSC{Float64,Int}
    bolt_law::SpringLaw
    link_law::SpringLaw               # contact stress – closure law
    link_area::Vector{Float64}
    bolt_dmax::Vector{Float64}
    link_dmax::Vector{Float64}
end

function Model(inp::BasePlateInput; mesh::Mesh = generate_mesh(inp))
    shell = ShellModel(mesh, inp)
    Ks = shell.K_lin + beam_stiffness(mesh, inp) + support_stiffness(mesh, inp)
    K_const = blockdiag(Ks, concrete_stiffness(mesh, inp))
    ncb = length(mesh.bp_nodes)
    area = inp.link_area == :average ? fill(sum(mesh.bp_area) / ncb, ncb) :
           inp.link_area == :tributary ? copy(mesh.bp_area) : throw(ArgumentError("link_area must be :average or :tributary"))
    return Model(inp, mesh, shell, K_const,
                 elastic_hardening_spring(inp.bolt_stiffness, inp.bolt_yield, inp.bolt_hardening),
                 elastic_hardening_spring(inp.contact_stiffness, inp.contact_yield, inp.contact_hardening),
                 area, zeros(length(mesh.bolt_nodes)), zeros(ncb))
end

bolt_dof(m::Mesh, i::Int) = 6 * (m.bolt_nodes[i] - 1) + 3
plate_dof(m::Mesh, i::Int) = 6 * (m.bp_nodes[i] - 1) + 3
concrete_dof(m::Mesh, i::Int) = 6 * nnodes(m) + i

"bolt force (tension positive) and tangent of bolt `i` for the displacements `u`"
bolt_response(model::Model, u, i; uncond = false) = uncond ? (model.bolt_law.k0 * u[bolt_dof(model.mesh, i)], model.bolt_law.k0) :
                                                    spring_response(model.bolt_law, u[bolt_dof(model.mesh, i)], model.bolt_dmax[i])

"contact force (compression positive) and tangent of link `i`"
function link_response(model::Model, u, i; uncond = false)
    δ = u[concrete_dof(model.mesh, i)] - u[plate_dof(model.mesh, i)]       # closure
    A = model.link_area[i]
    uncond && return model.link_law.k0 * A * δ, model.link_law.k0 * A
    s, k = spring_response(model.link_law, δ, model.link_dmax[i])
    return s * A, k * A
end

"""
    internal_forces(model, u; uncond = false) -> (f, K)

Internal force vector and tangent stiffness matrix for the displacements `u`. With `uncond` all
bolt springs and contact links act with their initial stiffness in both directions (used for the
very first iteration, where nothing is in contact yet).
"""
function internal_forces(model::Model, u::AbstractVector; uncond::Bool = false)
    m = model.mesh
    f = model.K_const * u
    assemble_nonlinear!(f, model.shell, u)
    I = Int[]; J = Int[]; V = Float64[]
    for i in eachindex(m.bolt_nodes)
        fb, kb = bolt_response(model, u, i; uncond)
        a = bolt_dof(m, i)
        f[a] += fb
        kb != 0 && (push!(I, a); push!(J, a); push!(V, kb))
    end
    for i in eachindex(m.bp_nodes)
        fc, kc = link_response(model, u, i; uncond)
        a = plate_dof(m, i); c = concrete_dof(m, i)
        f[a] -= fc; f[c] += fc
        if kc != 0
            append!(I, (a, c, a, c)); append!(J, (a, c, c, a)); append!(V, (kc, kc, -kc, -kc))
        end
    end
    n = length(u)
    K = model.K_const + model.shell.K_nl + sparse(I, J, V, n, n)
    return f, K
end

function factorize_tangent(K)
    try
        return cholesky(Symmetric(K))
    catch err
        err isa PosDefException || err isa ArgumentError || rethrow()
        return lu(K)
    end
end

"accept the current state as converged: element states and spring history"
function commit!(model::Model, u)
    commit!(model.shell)
    m = model.mesh
    for i in eachindex(m.bolt_nodes)
        model.bolt_dmax[i] = max(model.bolt_dmax[i], u[bolt_dof(m, i)])
    end
    for i in eachindex(m.bp_nodes)
        model.link_dmax[i] = max(model.link_dmax[i], u[concrete_dof(m, i)] - u[plate_dof(m, i)])
    end
    return model
end

"load vector with the force/moment `p = (Fx, Fy, Fz, Mx, My, Mz)` at the master node"
function master_load(m::Mesh, p)
    lv = zeros(ndofs(m))
    lv[nodedofs(m.master)] .= p
    return lv
end

"rotation of the upright between the two flange measurement points, and their mean vertical displacement"
function rotation(model::Model, u)
    m = model.mesh
    h = centerline(model.inp.geometry)[1]
    w1 = u[6*(m.rot_nodes[1]-1)+3]; w2 = u[6*(m.rot_nodes[2]-1)+3]
    return (w2 - w1) / h, (w1 + w2) / 2
end

"""
    Frame

State at the end of a converged increment: `step` (1 = axial, 2 = moment), axial force `N`
(compression positive), base moment `M`, rotation `θ` (absolute value, as in `bp_Barcelona.m`),
mean vertical displacement `w` of the two measurement points, Newton iterations, displacement
vector `u`, bolt forces, contact stresses at the base plate nodes and element results
(see [`element_fields`](@ref)).
"""
struct Frame
    step::Int
    N::Float64
    M::Float64
    θ::Float64
    w::Float64
    iterations::Int
    u::Vector{Float64}
    bolt_force::Vector{Float64}
    contact_stress::Vector{Float64}
    fields::Dict{Symbol,Vector{Float64}}
end

function Frame(model::Model, step, N, M, iters, u)
    m = model.mesh
    θ, w = rotation(model, u)
    bolt = [bolt_response(model, u, i)[1] for i in eachindex(m.bolt_nodes)]
    contact = [link_response(model, u, i)[1] / m.bp_area[i] for i in eachindex(m.bp_nodes)]
    return Frame(step, N, M, abs(θ), w, iters, copy(u), bolt, contact, element_fields(model, u))
end

"""
    Result

`model`, the converged `frames` (the first one is the unloaded state), and `status`: `:completed`,
`:max_rotation`, `:max_steps` or `:failed` (no convergence with the smallest increment).
"""
struct Result
    model::Model
    frames::Vector{Frame}
    status::Symbol
end

"""
    base_moment(model, frame; second_order = true)

Moment at the base plate. `frame.M` is the first-order moment of the lateral load (`V Ht`); with
`second_order` the moment of the axial force on the lateral displacement of the loaded
cross-section relative to the anchor bolts is added (`N δ`), as in the moment definition of the
Barcelona tests. Meaningful only with geometrically nonlinear upright elements.
"""
function base_moment(model::Model, f::Frame; second_order::Bool = true)
    second_order || return f.M
    m = model.mesh
    δ = f.u[6*(m.master-1)+2] - sum(f.u[6*(i-1)+2] for i in m.bolt_nodes) / length(m.bolt_nodes)
    return f.M + f.N * δ
end

"""
    moment_rotation(result; second_order = false) -> (θ, M)

Moment–rotation curve of step 2 (starting at the end of the axial preload); see [`base_moment`](@ref).
"""
function moment_rotation(r::Result; second_order::Bool = false)
    i0 = findlast(f -> f.step == 1, r.frames)
    fr = r.frames[i0:end]
    return [f.θ for f in fr], [base_moment(r.model, f; second_order) for f in fr]
end

# Newton iteration on the residual f_int(u) − fext(λ); in displacement control the load factor is
# an unknown and the constraint c·u = q closes the system (bordered solve with one factorization).
function newton!(model::Model, u0, fext0, pattern, λ0; λ = nothing, q = nothing, c = nothing,
                 rtol, atol, maxiter, uncond_first = false, verbose = false)
    u = copy(u0)
    λk = λ === nothing ? λ0 : λ
    restore!(model.shell)
    for it in 1:maxiter
        f, K = internal_forces(model, u; uncond = uncond_first && it == 1)
        r = f .- fext0 .- λk .* pattern
        rn = norm(r)
        fref = max(norm(fext0 .+ λk .* pattern), norm(f), atol)
        verbose && @printf("      it %2d  |r| = %.3e  λ = %.5g\n", it, rn, λk)
        if it > 1 && rn <= rtol * fref && (q === nothing || abs(dot(c, u) - q) <= 1e-9 * max(abs(q), 1e-12))
            return u, λk, it - 1, true
        end
        (isfinite(rn) && rn < 1e12 * fref) || break
        F = factorize_tangent(K)
        if q === nothing
            u .-= F \ r
        else
            dr = F \ r
            dl = F \ pattern
            dλ = (q - dot(c, u) + dot(c, dr)) / dot(c, dl)
            u .+= dλ .* dl .- dr
            λk += dλ
        end
    end
    return u, λk, maxiter, false
end

"""
    analyze(model; N, Mmax, θmax = Inf, kwargs...) -> Result
    analyze(inp::BasePlateInput; kwargs...)

Two-step analysis: the axial compression `N` (kip) is applied in load control, then the lateral
load pattern of `bp_Barcelona.m` (`V = M/Ht` at the master node at height `H` plus the moment
`M − V H`, bending about the upright's axis of symmetry) is increased until the base moment reaches
`Mmax` or the rotation reaches `θmax`.

Keywords: `axial_increments = 5`, `dM0 = Mmax/100` (first moment increment), `dM_max = Mmax/20`,
`dθ_max = θmax/40` (largest increments), `rtol = 1e-6` (residual force tolerance relative to the
force level), `maxiter = 20`, `max_steps = 400`, `verbose = true`.
"""
function analyze(model::Model; N::Real, Mmax::Real, θmax::Real = Inf, axial_increments::Int = 5,
                 dM0::Real = Mmax / 100, dM_max::Real = Mmax / 20, dθ_max::Real = θmax / 40,
                 rtol::Real = 1e-6, atol::Real = 1e-9, maxiter::Int = 20, max_steps::Int = 400, verbose::Bool = true)
    m = model.mesh
    H, Ht, _ = rounded_heights(model.inp)
    n = ndofs(m)
    u = zeros(n)
    frames = [Frame(model, 1, 0.0, 0.0, 0, u)]
    P1 = master_load(m, (0.0, 0.0, -float(N), 0.0, 0.0, 0.0))
    zero_load = zeros(n)

    # ---- step 1: axial force, load control ----------------------------------------------
    λ = 0.0; dλ = 1 / axial_increments; first = true
    while λ < 1 - 1e-12
        λt = min(λ + dλ, 1.0)
        ut, _, it, ok = newton!(model, u, zero_load, P1, λt; λ = λt, rtol, atol, maxiter, uncond_first = first)
        if !ok
            dλ /= 4
            dλ < 1e-6 && return Result(model, frames, :failed)
            continue
        end
        first = false
        u = ut; λ = λt
        commit!(model, u)
        push!(frames, Frame(model, 1, λ * N, 0.0, it, u))
        verbose && @printf("  axial   N = %7.3f  w = %10.3e  iterations = %d\n", λ * N, frames[end].w, it)
    end

    # ---- step 2: moment, displacement control on q = c·u -----------------------------------
    # unit base moment: V = 1/Ht in +y at height H, moment −(1 − H/Ht) about x
    P2 = master_load(m, (0.0, 1 / Ht, 0.0, -(1 - H / Ht), 0.0, 0.0))
    c = P2
    q = dot(c, u)
    _, K = internal_forces(model, u)
    restore!(model.shell)
    dq_dM = dot(c, factorize_tangent(K) \ P2)        # tangent compliance of the control displacement
    dq = dM0 * dq_dM
    M = 0.0; θ = frames[end].θ; status = :max_steps
    nfail = 0
    for _ in 1:max_steps
        ut, Mt, it, ok = newton!(model, u, P1, P2, M; q = q + dq, c, rtol, atol, maxiter)
        if !ok
            dq /= 4; nfail += 1
            if nfail > 8
                status = :failed
                break
            end
            continue
        end
        nfail = 0
        commit!(model, ut)
        fr = Frame(model, 2, float(N), Mt, it, ut)
        push!(frames, fr)
        verbose && @printf("  moment  M = %8.3f  θ = %9.6f  iterations = %d\n", Mt, fr.θ, it)
        dM = Mt - M; dθ = fr.θ - θ
        u = ut; q += dq; M = Mt; θ = fr.θ
        if M >= Mmax * (1 - 1e-9)
            status = :completed
            break
        elseif θ >= θmax
            status = :max_rotation
            break
        end
        # next increment: grow or shrink with the iteration count, limited in ΔM and Δθ
        grow = it <= 5 ? 1.5 : it <= 10 ? 1.2 : it <= 15 ? 1.0 : 0.5
        lim = Inf
        dM > 0 && (lim = min(lim, dq * min(dM_max, Mmax - M + 1e-6 * Mmax) / dM))
        (dθ > 0 && isfinite(dθ_max)) && (lim = min(lim, dq * dθ_max / dθ))
        dq = min(dq * grow, lim)
    end
    return Result(model, frames, status)
end

analyze(inp::BasePlateInput; mesh::Mesh = generate_mesh(inp), kwargs...) = analyze(Model(inp; mesh); kwargs...)

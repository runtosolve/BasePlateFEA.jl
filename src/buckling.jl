# --------------------------------------------------------------------------------------------
# Linearized buckling modes of the model, used as shapes of the initial geometric imperfection.
#
# The geometric stiffness comes from the nonlinear shell elements themselves: it is the derivative
# of their tangent stiffness along the elastic reference displacements (central difference), so it
# contains the initial stress and the initial displacement terms of the Total Lagrangian elements.
# Bolt springs and contact links act with their initial stiffness in both directions.
# --------------------------------------------------------------------------------------------

"""
    buckling_modes(model; N, M, nmodes = 4, subspace = 3nmodes, iterations = 60, tol = 1e-6) -> (λ, Φ)

Buckling load factors `λ` (ascending, positive) and modes `Φ` (columns, global dofs) for the
reference loading axial force `N` plus base moment `M` (loads at buckling are `λ N`, `λ M`).
Only the geometrically nonlinear shell elements contribute to the geometric stiffness, so the
upright steel must be plastic (or `nonlinear_elements = :all`). The model state is not changed.
"""
function buckling_modes(model::Model; N::Real, M::Real, nmodes::Int = 4, subspace::Int = 3nmodes,
                        iterations::Int = 60, tol::Real = 1e-6, fd::Real = 0.01)
    m = model.mesh; s = model.shell
    has_nonlinear(s) || throw(ArgumentError("buckling_modes needs geometrically nonlinear shell elements"))
    H, Ht, _ = rounded_heights(model.inp)
    n = ndofs(m)
    restore!(s)
    _, K0 = internal_forces(model, zeros(n); uncond = true)
    K0 = sparse(Symmetric(K0))
    F = cholesky(Symmetric(K0))
    P = master_load(m, (0.0, M / Ht, -float(N), -M * (1 - H / Ht), 0.0, 0.0))
    uref = F \ P
    function Knl(u)
        restore!(s)
        f = zeros(n)
        assemble_nonlinear!(f, s, u); assemble_nonlinear!(f, s, u)   # second pass updates the internal bubble dofs
        K = copy(s.K_nl)
        restore!(s)
        return K
    end
    Kg = (Knl(fd .* uref) - Knl(-fd .* uref)) ./ (2fd)
    Kg = (Kg + Kg') ./ 2

    # subspace iteration for the largest μ = 1/λ of  −Kg φ = μ K0 φ
    X = randn(MersenneTwisterSeed(1), n, subspace)
    μ = zeros(subspace); μold = fill(Inf, subspace)
    for it in 1:iterations
        Y = F \ (-(Kg * X))
        A = Symmetric(Y' * (K0 * Y)); B = Symmetric(-(Y' * (Kg * Y)))
        E = eigen(B, A)
        p = sortperm(E.values; rev = true)
        μ = E.values[p]
        X = Y * E.vectors[:, p]
        maximum(abs.(μ[1:nmodes] .- μold[1:nmodes]) ./ abs.(μ[1:nmodes])) < tol && break
        μold = copy(μ)
    end
    all(μ[1:nmodes] .> 0) || @warn "fewer than $nmodes positive buckling load factors were found"
    return 1 ./ μ[1:nmodes], X[:, 1:nmodes]
end

MersenneTwisterSeed(seed) = Random.MersenneTwister(seed)

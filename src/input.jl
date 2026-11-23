# --------------------------------------------------------------------------------------------
# Model input. The defaults reproduce the settings of S. Ádány's `bp_Barcelona.m`
# (to_Cris_Nov6_2026). Units are kip and inch throughout.
# --------------------------------------------------------------------------------------------

"""
    BasePlateGeometry(; W, L, t, Da, Dc, Sa, R, h, b, c, tu, hole_sectors = 16)

Lipped-C upright on a rectangular base plate with two anchor bolts, in the frame of
`bp_Barcelona.m`: the plate occupies `0 ≤ x ≤ L`, `0 ≤ y ≤ W`, `z = 0`; the upright's axis of
symmetry is the line `y = W/2`, its web is at `x = Dc` and its flanges run in `+x`; the two bolt
holes (radius `R`) are at `x = L − Da`, `y = W/2 ± Sa/2`.

`h`, `b`, `c` are the out-to-out web depth, flange width and lip length of the upright and `tu` its
thickness; the model uses the centerline dimensions `h − tu`, `b − tu`, `c − tu/2`, with the web
centerline at `x = Dc`.
"""
Base.@kwdef struct BasePlateGeometry
    W::Float64
    L::Float64
    t::Float64
    Da::Float64
    Dc::Float64
    Sa::Float64
    R::Float64
    h::Float64
    b::Float64
    c::Float64
    tu::Float64
    hole_sectors::Int = 16
end

"centerline dimensions `(h, b, c)` of the upright"
centerline(g::BasePlateGeometry) = (g.h - g.tu, g.b - g.tu, g.c - g.tu / 2)

"""
    barcelona_geometry(type)

Geometry of the Barcelona base plate types `:A`, `:B`, `:C` as defined in `bp_Barcelona.m`
(November 2026).
"""
function barcelona_geometry(type::Symbol)
    if type == :A
        R = 0.625 / 2
        return BasePlateGeometry(W = 8.0, L = 5.0, t = 0.375, Da = 2.5, Dc = 1.0, Sa = 8 - 2 * (0.625 + R), R = R,
                                 h = 3.0, b = 3.0, c = 0.75, tu = 0.1)
    elseif type == :B
        return BasePlateGeometry(W = 5.0, L = 7.0, t = 0.372, Da = 0.75, Dc = 0.25, Sa = 3.5, R = 0.375,
                                 h = 2 + 15 / 16, b = 3.0, c = 0.75, tu = 0.1016)
    elseif type == :C
        return BasePlateGeometry(W = 3.25, L = 5.0, t = 0.128, Da = 0.625, Dc = 0.125, Sa = 1.75, R = 0.3125,
                                 h = 3.0, b = 3.0, c = 0.75, tu = 0.1016)
    end
    throw(ArgumentError("unknown Barcelona base plate type $type (use :A, :B or :C)"))
end

"`Mmax0` of `bp_Barcelona.m`: the base moment at the end of the analysis is `Mmax0 + N` (kip-in, N in kip)"
barcelona_Mmax0(type::Symbol) = type == :A ? 40.0 : type == :B ? 20.0 : type == :C ? 7.0 :
                                throw(ArgumentError("unknown Barcelona base plate type $type"))

"""
    Steel(; E = 29500, ν = 0.3, fy = Inf, εp = [0.0], σ = [fy])

Shell steel. `fy = Inf` is linear elastic. Otherwise von Mises plasticity with the isotropic
hardening curve `(εp, σ)`: true stress `σ[i]` at equivalent plastic strain `εp[i]`, `εp[1] = 0`,
constant beyond the last point (see `IsotropicHardening` in Tri/QuadShellFiniteElement.jl).
"""
Base.@kwdef struct Steel
    E::Float64 = 29500.0
    ν::Float64 = 0.3
    fy::Float64 = Inf
    εp::Vector{Float64} = [0.0]
    σ::Vector{Float64} = [fy]
end

isplastic(s::Steel) = isfinite(s.fy)

"""
    bilinear_steel(fy, fu, εu; E = 29500, ν = 0.3)

Elastic – linear hardening steel from the engineering yield stress `fy`, tensile strength `fu` and
the engineering strain at `fu` (`εu`), converted to true stress and logarithmic plastic strain;
perfectly plastic beyond the tensile strength.
"""
function bilinear_steel(fy, fu, εu; E = 29500.0, ν = 0.3)
    σu = fu * (1 + εu)
    εpu = log(1 + εu) - σu / E
    (σu > fy && εpu > 0) || return Steel(; E, ν, fy, εp = [0.0], σ = [fy])
    return Steel(; E, ν, fy, εp = [0.0, εpu], σ = [fy, σu])
end

"""
    yun_gardner_steel(fy, fu; E = 29500, ν = 0.3, εsh, εu, npoints = 25)

Steel with a yield plateau and rounded strain hardening: the quad-linear-type model of Yun &
Gardner (2017) for hot-rolled steel, `StressStrainCurves.YunGardner2017`, converted to true stress
and logarithmic plastic strain; perfectly plastic beyond the tensile strength. By default the
strain at the tensile strength and the strain hardening strain are the predictive expressions of
that paper, `εu = 0.6 (1 − fy/fu) ≥ 0.06` and `εsh = 0.1 fy/fu − 0.055` within `[0.015, 0.03]`.
"""
function yun_gardner_steel(fy, fu; E = 29500.0, ν = 0.3, εu = max(0.6 * (1 - fy / fu), 0.06),
                           εsh = clamp(0.1 * fy / fu - 0.055, 0.015, 0.03), npoints::Int = 25)
    εy = fy / E
    # strain points: yield, end of the plateau, then graded towards the start of the hardening branch
    ε = vcat(εy, [εsh + (εu - εsh) * ((i - 1) / (npoints - 1))^2 for i in 1:npoints])
    ε[end] = εu
    curve = StressStrainCurves.YunGardner2017(E, fy, fu, εsh, εu, ε)
    σt = collect(curve.σ) .* (1 .+ ε)
    εp = log.(1 .+ ε) .- σt ./ E
    εp[1] = 0.0
    keep = [1]
    for i in 2:length(εp)                       # strictly increasing plastic strain, non-decreasing stress
        (εp[i] > εp[keep[end]] + 1e-9 && σt[i] >= σt[keep[end]]) && push!(keep, i)
    end
    return Steel(; E, ν, fy = σt[1], εp = εp[keep], σ = σt[keep])
end

"""
    SpringLaw(d, f)

Piecewise-linear force–displacement law of a one-sided spring as in `MAT_SPR_table`: the points
`(d[i], f[i])` with `d[1] < 0`, `f[1] = 0`, `(d[2], f[2]) = (0, 0)`; no force for `d < 0`, constant
force beyond the last point. Unloading from the largest displacement reached so far follows the
initial stiffness down to zero force (`spr_update.m`).
"""
struct SpringLaw
    d::Vector{Float64}
    f::Vector{Float64}
    k0::Float64
    function SpringLaw(d, f)
        length(d) == length(f) ≥ 3 || throw(ArgumentError("a spring law needs at least 3 points"))
        return new(collect(Float64, d), collect(Float64, f), (f[3] - f[2]) / (d[3] - d[2]))
    end
end

"elastic (stiffness `k`) – yield (`fy`) – hardening spring: the 4-point law of `bp_Barcelona.m`, hardening to `fmult × fy` at `101 dy`"
function elastic_hardening_spring(k, fy, fmult)
    dy = fy / k
    return SpringLaw([-1000.0, 0.0, dy, 101 * dy], [0.0, 0.0, fy, fmult * fy])
end

"""
    BasePlateInput(; geometry, kwargs...)

Complete model definition. All keyword defaults are the values of `bp_Barcelona.m` (November 2026)
except where noted.

Mesh and upright
- `fe_size = 0.15`: target element size
- `H = 3.0`: modelled upright height (rounded to a multiple of `fe_size`)
- `Ht = 24.0`: height of the lateral load above the base plate in the test (rounded likewise);
  the load at the top of the model is `V = M/Ht` with the moment `M − V H`
- `H_fine = Inf`, `fe_size_upper = fe_size`: the upright walls are meshed with `fe_size` up to the
  height `H_fine` and with the element height `fe_size_upper` above it (for tall uprights)
- `upright_plastic_height = Inf`: upright elements above this height keep an elastic steel (they
  are still geometrically nonlinear if the upright steel is plastic)
- `top_restraint = :free`: `:planar` holds the loaded cross-section in the plane of bending (no
  out-of-plane translation, no twist, no rotation about the bending-plane axis), as a test rig
  with planar hinges does; needed for tall uprights, whose open section is otherwise free to twist
- `Hm = H`: height at which the rotation is measured, `ym_ratio = 0.5`: position of the two
  measurement points along the flanges, as a fraction of the flange width from the web

Materials
- `plate_steel`, `upright_steel :: Steel` (elastic by default)
- `Ec = 3410`, `νc = 0.2`: concrete half-space
- `bolt_stiffness = 120` (kip/in), `bolt_yield = 3500` (kip), hardening to `2 fy` at `101 dy`
- `contact_stiffness = 5000` (ksi/in), `contact_yield = 50` (ksi), hardening to `11 fy`
- `friction_stiffness = 0.1`, `hole_stiffness = 1.0` (kip/in per node, in-plane x and y)
- `top_mpc = (A = 100, Iy = 1e4, Iz = 1e4, J = 1e2)`, `hole_mpc = (A = 1, Iy = 1, Iz = 1, J = 1)`:
  section properties of the rigid-link beams (steel `E`, `G`)
- `hole_thickness_factor = 0.01`: thickness multiplier of the plate elements inside the holes

Shell element options (Tri/QuadShellFiniteElement.jl)
- `Cs_tri = 0.0`, `Cs_quad = 0.0`, `drilling = :penalty`: Ádány's MATLAB elements. The packages'
  own defaults are `Cs_tri = 0.2`, `Cs_quad = 0.1`, `drilling = :hughes_brezzi`.
- `nlayers = 5`, `layer_rule = :simpson`: through-thickness integration of plastic sections
- `nonlinear_elements = :auto`: `:auto` uses the geometrically nonlinear elastic–plastic elements
  only for parts with a plastic steel, `:all` for every shell element, `:none` never (linear)

Concrete discretization
- `patch_factor = 0.1`: side of the square Boussinesq patch as a fraction of √(nodal area)
- `concrete_threshold = 0.01`: entries of the concrete stiffness matrix below this fraction of
  the largest entry are dropped
- `link_area = :average`: nodal area used to turn the contact law into forces (`:average` as in
  the MATLAB code, or `:tributary`)
"""
Base.@kwdef struct BasePlateInput
    geometry::BasePlateGeometry
    fe_size::Float64 = 0.15
    H::Float64 = 3.0
    Ht::Float64 = 24.0
    H_fine::Float64 = Inf
    fe_size_upper::Float64 = fe_size
    upright_plastic_height::Float64 = Inf
    top_restraint::Symbol = :free
    Hm::Float64 = H
    ym_ratio::Float64 = 0.5
    plate_steel::Steel = Steel()
    upright_steel::Steel = Steel()
    Ec::Float64 = 3410.0
    νc::Float64 = 0.2
    bolt_stiffness::Float64 = 120.0
    bolt_yield::Float64 = 3500.0
    bolt_hardening::Float64 = 2.0
    contact_stiffness::Float64 = 5000.0
    contact_yield::Float64 = 50.0
    contact_hardening::Float64 = 11.0
    friction_stiffness::Float64 = 0.1
    hole_stiffness::Float64 = 1.0
    top_mpc::NamedTuple{(:A, :Iy, :Iz, :J),NTuple{4,Float64}} = (A = 100.0, Iy = 1e4, Iz = 1e4, J = 1e2)
    hole_mpc::NamedTuple{(:A, :Iy, :Iz, :J),NTuple{4,Float64}} = (A = 1.0, Iy = 1.0, Iz = 1.0, J = 1.0)
    hole_thickness_factor::Float64 = 0.01
    Cs_tri::Float64 = 0.0
    Cs_quad::Float64 = 0.0
    drilling::Symbol = :penalty
    drilling_gamma::Float64 = 1.0
    nlayers::Int = 5
    layer_rule::Symbol = :simpson
    nonlinear_elements::Symbol = :auto
    patch_factor::Float64 = 0.1
    concrete_threshold::Float64 = 0.01
    link_area::Symbol = :average
end

"modelled upright height, lever arm and measurement height rounded to the mesh as in `bp_Barcelona.m`"
function rounded_heights(inp::BasePlateInput)
    a = inp.fe_size
    H = round(inp.H / a) * a
    Ht = round(inp.Ht / a) * a
    Hm = min(round(inp.Hm / a) * a, H)
    return H, Ht, Hm
end

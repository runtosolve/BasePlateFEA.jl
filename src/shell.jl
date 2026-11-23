# --------------------------------------------------------------------------------------------
# Shell elements. The 3-node triangle and the 4-node quadrilateral of S. Ádány's MATLAB code are
# the elements of TriShellFiniteElement.jl and QuadShellFiniteElement.jl (with `Cs = 0` and
# `drilling = :penalty`), so the element matrices come from those packages:
#
#  * linear elements: `local_elastic_stiffness_matrix!`, assembled once;
#  * nonlinear elements: the geometrically nonlinear, layered J2-plasticity elements of the two
#    packages (`tri_internal_force!`, `quad_internal_force!`), evaluated at every iteration. The
#    element data are built here element by element, as in the packages' `PlasticShellModel`
#    constructors, because the base plate model mixes triangles and quadrilaterals of different
#    thickness and material.
#
# Element dofs are in nodal block order [ux uy uz θx θy θz] per node, in global axes.
# --------------------------------------------------------------------------------------------

const TRI_QR1 = QuadratureRule{RefTriangle}(1)
const TRI_QR3 = QuadratureRule{RefTriangle}(2)
const QUAD_QRM = QuadratureRule{RefQuadrilateral}(2)
const QUAD_QRB = QuadratureRule{RefQuadrilateral}(3)

element_coordinates(m::Mesh, e::Int) = [Vec{3}((m.x[i, 1], m.x[i, 2], m.x[i, 3])) for i in m.conn[e]]
element_dofs(m::Mesh, e::Int) = reduce(vcat, [collect(nodedofs(i)) for i in m.conn[e]])
part_steel(inp::BasePlateInput, part::Int) = part == 1 ? inp.plate_steel : inp.upright_steel

"global elastic stiffness matrix (18×18 or 24×24) of a linear shell element"
function elastic_element_stiffness(xg, E, ν, t, inp::BasePlateInput)
    if length(xg) == 3
        T = TSF.calculation_rotation_matrix(xg)
        xl = TSF.global_nodal_coords_to_planar_coords(xg, T)
        k = TSF.local_elastic_stiffness_matrix!(TRI_QR1, TRI_QR3, TSF.IP3(), TSF.IP6(), E, ν, t, xl;
                                                Cs = inp.Cs_tri, drilling = inp.drilling, drilling_gamma = inp.drilling_gamma)
        Te = TSF.rotation_matrix_for_element_stiffness_drilling(T)
    else
        T = QSF.calculation_rotation_matrix(xg)
        xl = QSF.global_nodal_coords_to_planar_coords(xg, T)
        k = QSF.local_elastic_stiffness_matrix!(QUAD_QRM, QUAD_QRB, QSF.IP4(), QSF.IP6(), E, ν, t, xl;
                                                Cs = inp.Cs_quad, drilling = inp.drilling, drilling_gamma = inp.drilling_gamma)
        Te = QSF.rotation_matrix_for_element_stiffness_drilling(T)
    end
    return Te * k * Te'
end

# element data of the nonlinear triangle, as in TSF.PlasticShellModel
function tri_element_data(xg, sec::TSF.ShellSection, inp::BasePlateInput)
    E = sec.mat.E; ν = sec.mat.ν; t = sec.t
    ip3 = TSF.IP3(); ip6 = TSF.IP6()
    cv = CellValues(TRI_QR1, ip3, ip3)
    T = TSF.calculation_rotation_matrix(xg)
    xl = TSF.global_nodal_coords_to_planar_coords(xg, T)
    reinit!(cv, xl)
    dNdx = [shape_gradient(cv, 1, i) for i in 1:3]
    area = getdetJdV(cv, 1)
    ke_m, ke_b, ke_s = TSF.local_elastic_pieces(TRI_QR1, TRI_QR3, ip3, ip6, E, ν, t, xl; Cs = inp.Cs_tri)
    induv = [1 2 6 7 11 12]; indwt = [3 4 5 8 9 10 13 14 15]
    ind18 = [1:5; 7:11; 13:17]
    k15 = zeros(15, 15)
    k15[induv, induv] = ke_m
    k15[indwt, indwt] = ke_s + ke_b
    k15_s = zeros(15, 15)
    k15_s[indwt, indwt] = ke_s
    k_el = zeros(18, 18)
    k_el[ind18, ind18] = k15_s
    TSF.add_drilling_stiffness!(k_el, k15, TRI_QR3, ip3, xl, E, ν, t, inp.drilling, inp.drilling_gamma)
    Te = TSF.rotation_matrix_for_element_stiffness_drilling(T)
    nl = TSF.nlayers(sec)
    return TSF.TriShellElementData(Matrix(T), Te, dNdx, area, k_el, [TSF.LayerState() for _ in 1:nl], [TSF.LayerState() for _ in 1:nl])
end

# element data of the nonlinear quadrilateral, as in QSF.PlasticShellModel
function quad_element_data(xg, sec::QSF.ShellSection, inp::BasePlateInput, cv_m, cv_b)
    E = sec.mat.E; ν = sec.mat.ν; t = sec.t
    nq = getnquadpoints(cv_b)
    T = QSF.calculation_rotation_matrix(xg)
    xl = QSF.global_nodal_coords_to_planar_coords(xg, T)
    ke_m, ke_b, ke_s = QSF.local_elastic_pieces(cv_m, cv_b, E, ν, t, xl; Cs = inp.Cs_quad)
    dNdx = [shape_gradient(cv_b, q, i) for i in 1:6, q in 1:nq]
    dV = [getdetJdV(cv_b, q) for q in 1:nq]
    ks = zeros(30, 30)
    iwt = [5 * (i - 1) + c for i in 1:6 for c in 3:5]
    ks[iwt, iwt] = ke_s
    k20 = zeros(20, 20)
    k20[QSF.IND_UV_20, QSF.IND_UV_20] = QSF.static_condensation(ke_m, 2 * QSF.NNODES)
    k20[QSF.IND_WT_20, QSF.IND_WT_20] = QSF.static_condensation(ke_b + ke_s, 3 * QSF.NNODES)
    kd = zeros(24, 24)
    QSF.add_drilling_stiffness!(kd, cv_m, k20, E, ν, t, inp.drilling, inp.drilling_gamma)
    Te = QSF.rotation_matrix_for_element_stiffness_drilling(T)
    nl = QSF.nlayers(sec)
    states = [QSF.LayerState() for _ in 1:nl, _ in 1:nq]
    return QSF.QuadShellElementData(Matrix(T), Te, dNdx, dV, ks, kd, states, copy(states),
                                    zeros(10), zeros(24), zeros(10), zeros(10, 24), zeros(10, 10),
                                    zeros(10), zeros(24), zeros(10), zeros(10, 24), zeros(10, 10))
end

"""
    ShellModel

Shell part of the model: the assembled stiffness of the linear elements (`K_lin`, steel dofs) and
the data of the nonlinear elements. `tri`/`quad` hold the mesh element numbers of the nonlinear
triangles/quadrilaterals, `*_data` their element data and `*_sec` their sections.
"""
struct ShellModel
    K_lin::SparseMatrixCSC{Float64,Int}
    tri::Vector{Int}
    tri_data::Vector{TSF.TriShellElementData}
    tri_sec::Vector{TSF.ShellSection}
    quad::Vector{Int}
    quad_data::Vector{QSF.QuadShellElementData}
    quad_sec::Vector{QSF.ShellSection}
    ws::QSF.QuadWorkspace
    # assembly of the nonlinear elements into a matrix with a fixed sparsity pattern (all dofs)
    K_nl::SparseMatrixCSC{Float64,Int}
    pos::Vector{Int}                 # nzval position of every element matrix entry, element by element
    dofs::Vector{Vector{Int}}        # global dofs of the nonlinear elements (triangles first)
end

has_nonlinear(s::ShellModel) = !isempty(s.tri) || !isempty(s.quad)

function is_nonlinear_element(inp::BasePlateInput, m::Mesh, e::Int)
    inp.nonlinear_elements == :none && return false
    m.hole[e] && return false                      # the soft elements inside the holes stay linear elastic
    inp.nonlinear_elements == :all && return true
    inp.nonlinear_elements == :auto || throw(ArgumentError("nonlinear_elements must be :auto, :all or :none"))
    return isplastic(part_steel(inp, m.part[e]))
end

function ShellModel(m::Mesh, inp::BasePlateInput)
    ns = ndofs_steel(m); nt = ndofs(m)
    I = Int[]; J = Int[]; V = Float64[]
    tri = Int[]; quad = Int[]
    tri_data = TSF.TriShellElementData[]; quad_data = QSF.QuadShellElementData[]
    tri_sec = TSF.ShellSection[]; quad_sec = QSF.ShellSection[]
    tsecs = Dict{Tuple{Int,Float64},TSF.ShellSection}(); qsecs = Dict{Tuple{Int,Float64},QSF.ShellSection}()
    cv_m = CellValues(QUAD_QRM, QSF.IP6(), QSF.IP4()); cv_b = CellValues(QUAD_QRB, QSF.IP6(), QSF.IP4())
    for e in 1:nelements(m)
        xg = element_coordinates(m, e)
        st = part_steel(inp, m.part[e]); t = m.thickness[e]
        if is_nonlinear_element(inp, m, e)
            key = (m.part[e] == 1 ? 1 : 2, t)
            if m.part[e] > 1 && sum(q[3] for q in xg) / length(xg) > inp.upright_plastic_height
                st = Steel(; E = st.E, ν = st.ν, fy = 1e9, εp = [0.0], σ = [1e9])     # elastic, geometrically nonlinear
                key = (3, t)
            end
            if length(xg) == 3
                sec = get!(tsecs, key) do
                    TSF.ShellSection(TSF.J2PlaneStress(st.E, st.ν, TSF.IsotropicHardening(st.εp, st.σ)), t;
                                     nlayers = inp.nlayers, rule = inp.layer_rule)
                end
                push!(tri, e); push!(tri_sec, sec); push!(tri_data, tri_element_data(xg, sec, inp))
            else
                sec = get!(qsecs, key) do
                    QSF.ShellSection(QSF.J2PlaneStress(st.E, st.ν, QSF.IsotropicHardening(st.εp, st.σ)), t;
                                     nlayers = inp.nlayers, rule = inp.layer_rule)
                end
                push!(quad, e); push!(quad_sec, sec); push!(quad_data, quad_element_data(xg, sec, inp, cv_m, cv_b))
            end
        else
            ke = elastic_element_stiffness(xg, st.E, st.ν, t, inp)
            dofs = element_dofs(m, e)
            for b in eachindex(dofs), a in eachindex(dofs)
                push!(I, dofs[a]); push!(J, dofs[b]); push!(V, ke[a, b])
            end
        end
    end
    K_lin = sparse(I, J, V, ns, ns)

    # fixed sparsity pattern of the nonlinear elements
    dofs = [element_dofs(m, e) for e in vcat(tri, quad)]
    In = Int[]; Jn = Int[]
    for d in dofs, b in eachindex(d), a in eachindex(d)
        push!(In, d[a]); push!(Jn, d[b])
    end
    K_nl = sparse(In, Jn, zeros(length(In)), nt, nt)
    pos = Vector{Int}(undef, length(In))
    rows = rowvals(K_nl)
    for k in eachindex(In)
        r = nzrange(K_nl, Jn[k])
        pos[k] = first(r) + searchsortedfirst(view(rows, r), In[k]) - 1
    end
    return ShellModel(K_lin, tri, tri_data, tri_sec, quad, quad_data, quad_sec, QSF.QuadWorkspace(), K_nl, pos, dofs)
end

"""
    assemble_nonlinear!(f, s::ShellModel, u)

Add the internal forces of the nonlinear shell elements for the global displacements `u` to `f`
and fill `s.K_nl` with their tangent stiffness. The layer states of this evaluation are pending
until [`commit!`](@ref).
"""
function assemble_nonlinear!(f::AbstractVector, s::ShellModel, u::AbstractVector)
    nz = nonzeros(s.K_nl)
    fill!(nz, 0.0)
    k = 0
    ug = zeros(18); ul = zeros(18); fe = zeros(18); ke = zeros(18, 18); fg = zeros(18); kg = zeros(18, 18); tmp = zeros(18, 18)
    for (ie, el) in enumerate(s.tri_data)
        d = s.dofs[ie]
        @inbounds for a in 1:18
            ug[a] = u[d[a]]
        end
        mul!(ul, el.Te', ug)
        TSF.tri_internal_force!(fe, ke, el, s.tri_sec[ie], ul)
        mul!(fg, el.Te, fe)
        mul!(tmp, ke, el.Te'); mul!(kg, el.Te, tmp)
        @inbounds for a in 1:18
            f[d[a]] += fg[a]
        end
        @inbounds for b in 1:18, a in 1:18
            k += 1
            nz[s.pos[k]] += kg[a, b]
        end
    end
    ws = s.ws
    ug = ws.ug; ul = ws.ul; fe = ws.f24; ke = ws.k24; fg = ws.fg; kg = ws.kg; tmp = ws.tmp24
    nt = length(s.tri_data)
    for (ie, el) in enumerate(s.quad_data)
        d = s.dofs[nt+ie]
        @inbounds for a in 1:24
            ug[a] = u[d[a]]
        end
        mul!(ul, el.Te', ug)
        QSF.quad_internal_force!(fe, ke, el, s.quad_sec[ie], ws, ul)
        mul!(fg, el.Te, fe)
        mul!(tmp, ke, el.Te'); mul!(kg, el.Te, tmp)
        @inbounds for a in 1:24
            f[d[a]] += fg[a]
        end
        @inbounds for b in 1:24, a in 1:24
            k += 1
            nz[s.pos[k]] += kg[a, b]
        end
    end
    return f
end

"accept the pending states of the nonlinear elements as converged"
function commit!(s::ShellModel)
    for el in s.tri_data
        copyto!(el.states, el.pending)
    end
    for el in s.quad_data
        copyto!(el.states, el.pending)
        copyto!(el.α_c, el.α); copyto!(el.d_last_c, el.d_last); copyto!(el.rα_c, el.rα)
        copyto!(el.Kαd_c, el.Kαd); copyto!(el.Kαα_inv_c, el.Kαα_inv)
    end
    return s
end

"discard the pending states (after a rejected step)"
function restore!(s::ShellModel)
    for el in s.tri_data
        copyto!(el.pending, el.states)
    end
    for el in s.quad_data
        copyto!(el.pending, el.states)
        copyto!(el.α, el.α_c); copyto!(el.d_last, el.d_last_c); copyto!(el.rα, el.rα_c)
        copyto!(el.Kαd, el.Kαd_c); copyto!(el.Kαα_inv, el.Kαα_inv_c)
    end
    return s
end

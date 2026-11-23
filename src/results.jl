# --------------------------------------------------------------------------------------------
# Element results for output and visualization (main_stressout.m): stresses at the bottom, middle
# and top of the shell thickness and the equivalent plastic strain, one value per element.
# --------------------------------------------------------------------------------------------

_vm(s) = sqrt(max(s[1]^2 - s[1] * s[2] + s[2]^2 + 3 * s[3]^2, 0.0))

# gradients of the linear (triangle) or bilinear (quadrilateral, at the centroid) shape functions
# in the element's local frame; returns the frame, the gradients and the local nodal coordinates
function _centroid_gradients(xg)
    n = length(xg)
    T = n == 3 ? Matrix(TSF.calculation_rotation_matrix(xg)) : Matrix(QSF.calculation_rotation_matrix(xg))
    x0 = sum(xg) / n
    xl = [T' * collect(x - x0) for x in xg]
    if n == 3
        A2 = (xl[2][1] - xl[1][1]) * (xl[3][2] - xl[1][2]) - (xl[2][2] - xl[1][2]) * (xl[3][1] - xl[1][1])
        gx = [xl[2][2] - xl[3][2], xl[3][2] - xl[1][2], xl[1][2] - xl[2][2]] ./ A2
        gy = [xl[3][1] - xl[2][1], xl[1][1] - xl[3][1], xl[2][1] - xl[1][1]] ./ A2
    else
        dξ = [-1.0, 1.0, 1.0, -1.0] ./ 4; dη = [-1.0, -1.0, 1.0, 1.0] ./ 4
        J = [sum(dξ[i] * xl[i][1] for i in 1:4) sum(dξ[i] * xl[i][2] for i in 1:4)
             sum(dη[i] * xl[i][1] for i in 1:4) sum(dη[i] * xl[i][2] for i in 1:4)]
        Ji = inv(J)
        gx = [Ji[1, 1] * dξ[i] + Ji[1, 2] * dη[i] for i in 1:4]
        gy = [Ji[2, 1] * dξ[i] + Ji[2, 2] * dη[i] for i in 1:4]
    end
    return T, gx, gy
end

# rotate a plane stress vector from the element frame T to the plate frame Fp
function _to_plate_frame(s, T, Fp)
    Tse = Fp' * T
    c = Tse[1, 1]; sn = Tse[2, 1]
    return (c^2 * s[1] + sn^2 * s[2] - 2 * sn * c * s[3],
            sn^2 * s[1] + c^2 * s[2] + 2 * sn * c * s[3],
            sn * c * (s[1] - s[2]) + (c^2 - sn^2) * s[3])
end

"""
    element_fields(model, u) -> Dict{Symbol,Vector{Float64}}

One value per shell element (NaN for the soft elements inside the bolt holes):
`:vm_top`, `:vm_mid`, `:vm_bot` von Mises stress at `z = +t/2, 0, −t/2` (z along the element
normal: up for the base plate), `:sx`, `:sy`, `:txy` mid-surface (membrane) stresses in the plate
axes (base plate: global x, y; upright walls: along the wall, vertical), `:sx_top`, `:sy_top`,
`:sx_bot`, `:sy_bot` surface stresses in the same axes, and `:peeq` the largest equivalent plastic
strain through the thickness. Nonlinear elements report their converged layer stresses (averaged
over the integration points of a quadrilateral); linear elements are evaluated at the centroid.
"""
function element_fields(model::Model, u::AbstractVector)
    m = model.mesh; s = model.shell; inp = model.inp
    ne = nelements(m)
    names = (:vm_top, :vm_mid, :vm_bot, :sx, :sy, :txy, :sx_top, :sy_top, :sx_bot, :sy_bot, :peeq)
    out = Dict(k => fill(NaN, ne) for k in names)
    function store!(e, T, bot, mid, top, peeq)
        Fp = m.frames[m.part[e]]
        out[:vm_bot][e] = _vm(bot); out[:vm_mid][e] = _vm(mid); out[:vm_top][e] = _vm(top)
        pm = _to_plate_frame(mid, T, Fp); pt = _to_plate_frame(top, T, Fp); pb = _to_plate_frame(bot, T, Fp)
        out[:sx][e] = pm[1]; out[:sy][e] = pm[2]; out[:txy][e] = pm[3]
        out[:sx_top][e] = pt[1]; out[:sy_top][e] = pt[2]; out[:sx_bot][e] = pb[1]; out[:sy_bot][e] = pb[2]
        out[:peeq][e] = peeq
    end
    nonlinear = falses(ne)

    # nonlinear triangles
    ug = zeros(18); d = zeros(18)
    for (ie, el) in enumerate(s.tri_data)
        e = s.tri[ie]; nonlinear[e] = true
        sec = s.tri_sec[ie]
        ug .= u[s.dofs[ie]]
        mul!(d, el.Te', ug)
        dN = el.dNdx
        ux = uy = vx = vy = wx = wy = 0.0; κx = κy = κxy = 0.0
        for i in 1:3
            o = 6 * (i - 1); gx = dN[i][1]; gy = dN[i][2]
            ux += gx * d[o+1]; uy += gy * d[o+1]; vx += gx * d[o+2]; vy += gy * d[o+2]
            wx += gx * d[o+3]; wy += gy * d[o+3]
            κx += gx * d[o+5]; κy -= gy * d[o+4]; κxy += gy * d[o+5] - gx * d[o+4]
        end
        εm = SVector(ux + 0.5 * (ux^2 + vx^2 + wx^2), vy + 0.5 * (uy^2 + vy^2 + wy^2), uy + vx + ux * uy + vx * vy + wx * wy)
        σ = TSF.layer_stresses(sec, εm, SVector(κx, κy, κxy), el.states)
        nl = length(σ)
        store!(e, el.T, σ[1], σ[(nl+1)÷2], σ[nl], maximum(st.ε̄p for st in el.states))
    end

    # nonlinear quadrilaterals
    ug = zeros(24); d = zeros(24); U = zeros(30)
    nt = length(s.tri_data)
    for (ie, el) in enumerate(s.quad_data)
        e = s.quad[ie]; nonlinear[e] = true
        sec = s.quad_sec[ie]
        ug .= u[s.dofs[nt+ie]]
        mul!(d, el.Te', ug)
        for k in 1:20
            U[k] = d[QSF.MAP_20_TO_24[k]]
        end
        U[21:30] .= el.α_c
        nq = length(el.dV); nl = QSF.nlayers(sec)
        acc = [zero(SVector{3,Float64}) for _ in 1:nl]
        for q in 1:nq
            ux = uy = vx = vy = wx = wy = 0.0; κx = κy = κxy = 0.0
            for i in 1:6
                o = 5 * (i - 1); gx = el.dNdx[i, q][1]; gy = el.dNdx[i, q][2]
                ux += gx * U[o+1]; uy += gy * U[o+1]; vx += gx * U[o+2]; vy += gy * U[o+2]
                wx += gx * U[o+3]; wy += gy * U[o+3]
                κx += gx * U[o+5]; κy -= gy * U[o+4]; κxy += gy * U[o+5] - gx * U[o+4]
            end
            εm = SVector(ux + 0.5 * (ux^2 + vx^2 + wx^2), vy + 0.5 * (uy^2 + vy^2 + wy^2), uy + vx + ux * uy + vx * vy + wx * wy)
            σ = QSF.layer_stresses(sec, εm, SVector(κx, κy, κxy), view(el.states, :, q))
            for l in 1:nl
                acc[l] += σ[l] / nq
            end
        end
        store!(e, el.T, acc[1], acc[(nl+1)÷2], acc[nl], maximum(st.ε̄p for st in el.states))
    end

    # linear elements: centroid strains from the nodal displacements
    for e in 1:ne
        (nonlinear[e] || m.hole[e]) && continue
        xg = element_coordinates(m, e)
        T, gx, gy = _centroid_gradients(xg)
        st = part_steel(inp, m.part[e]); t = m.thickness[e]
        C = st.E / (1 - st.ν^2) * @SMatrix [1.0 st.ν 0.0; st.ν 1.0 0.0; 0.0 0.0 (1-st.ν)/2]
        εx = εy = γ = κx = κy = κxy = 0.0
        for (i, nd) in enumerate(m.conn[e])
            dofs = nodedofs(nd)
            tr = T' * u[dofs[1:3]]; ro = T' * u[dofs[4:6]]
            εx += gx[i] * tr[1]; εy += gy[i] * tr[2]; γ += gy[i] * tr[1] + gx[i] * tr[2]
            κx += gx[i] * ro[2]; κy -= gy[i] * ro[1]; κxy += gy[i] * ro[2] - gx[i] * ro[1]
        end
        εm = SVector(εx, εy, γ); κ = SVector(κx, κy, κxy)
        store!(e, T, C * (εm - t / 2 * κ), C * εm, C * (εm + t / 2 * κ), 0.0)
    end
    return out
end

"average the element values `v` to the nodes (ignoring NaN), for smooth contour plots; nodes without a value get NaN"
function nodal_average(m::Mesh, v::AbstractVector; parts = nothing)
    s = zeros(nnodes(m)); n = zeros(Int, nnodes(m))
    for (e, c) in enumerate(m.conn)
        (isnan(v[e]) || (parts !== nothing && !(m.part[e] in parts))) && continue
        for i in c
            s[i] += v[e]; n[i] += 1
        end
    end
    return [n[i] > 0 ? s[i] / n[i] : NaN for i in eachindex(s)]
end

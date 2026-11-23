module BasePlateFEA

"""
BasePlateFEA.jl - Column Base Plate Finite Element Analysis Package
Converted from MATLAB (S. Adany, latest update November 2026).

Models thin-walled cold-formed steel upright column connected to a steel base plate
resting on a concrete foundation, featuring:
  - 3-node triangular and 4-node quadrilateral shell elements (Mindlin plate theory)
  - DistMesh triangular meshing of the base plate around bolt holes & upright footprint
  - Multi-Point Constraint (MPC) 3D Euler-Bernoulli rigid beam links
  - Semi-infinite elastic half-space (Boussinesq) concrete foundation flexibility/stiffness
  - Compression-only nonlinear conditional contact links between plate and concrete
  - Nonlinear anchor bolt springs with piecewise constitutive response and unloading
  - Two-step incremental-iterative solver (axial preload + lateral/moment loading)
    with adaptive step sizing, stepback recovery, and Newton-Raphson iteration.
"""

using LinearAlgebra, SparseArrays, Statistics, Random
using GeometryBasics: Point3f, TriangleFace, Mesh

export Tolerances, FEData, HoleData
export cs_def, cs_prop, cs_data
export GL3, GL4, FE_props
export rotate3d, ct_3node_g2e, ct_4node_g2e, plate_transform
export collect_anchor, mesh4_rect, edges_by_plate
export delaunay_2d, deduplicate_points, distmesh
export HOLE_TABLE_entry, elem_divide_by_plate, four2three_by_plate
export main_modelgen, main_modelgen_dm
export ke_beam_EB, stiffmat_e_beam, stiffmat_concr
export dir_convert, collect_pl_data, collect_control, collect_bp, collect_edge, calc_bp_area, collect_spring
export spr_stif, spr_force, spr_update, spring_support, stifmat_cond, combine_Ke, force_back
export extend_lv, load_vec_point_MPC, lv_el3_edge, lv_el4_edge, load_vec_point, load_vec_UDL, load_vec_edge
export ke_uv_3n_nocondens_6dof, ke_uv_3n_condens_from_8to6dof, ke_wt_3n_condens_from_12to9dof_MIN
export ke_uv_4n_nocondens_8dof_num, ke_uv_4n_condens_from_12to8dof_num, ke_wt_4n_condens_from_18to12dof_num
export add_drill, stiffmat_e, fem_support, cond_support
export stress_uv_3n, stress_uv_4n, main_stressout
export FEanal_incr, FEanal_incr_old

# ─────────────────────────────────────────────────────────────────────────────
# TYPES
# ─────────────────────────────────────────────────────────────────────────────

mutable struct Tolerances
    zero_thickn::Float64
    zero_spring::Float64
    zero_dist::Float64
    zeroFElength::Float64
end

mutable struct FEData
    shape::Int
    nodenr::Int
    m_dof_ini::Int
    m_dof_final::Int
    m_ifcondens::Bool
    m_numint::Bool
    m_nGpe::Int
    b_theory::String
    b_dof_ini::Int
    b_dof_final::Int
    b_ifcondens::Bool
    b_numint::Bool
    b_nGpe::Int
    nonzer_ke::Int
    nonzer_lv::Int
    nGpg::Int
end

mutable struct HoleData
    coord::Vector{Vector{Float64}}
    radius::Vector{Float64}
    sector::Vector{Int}
    holes::Vector{Matrix{Float64}}
end

# ─────────────────────────────────────────────────────────────────────────────
# CROSS-SECTION
# ─────────────────────────────────────────────────────────────────────────────

function cs_def(h, b, c, tb, th, tc)
    cs_node = [
        1.0  b      -h/2+c
        2.0  b      -h/2
        3.0  0.0    -h/2
        4.0  0.0     h/2
        5.0  b       h/2
        6.0  b       h/2-c
    ]
    cs_elem = [
        1.0  1.0  2.0  tc
        2.0  2.0  3.0  tb
        3.0  3.0  4.0  th
        4.0  4.0  5.0  tb
        5.0  5.0  6.0  tc
    ]
    return cs_node, cs_elem
end

function cs_prop(cs_node::AbstractMatrix, cs_elem::AbstractMatrix)
    nelem = size(cs_elem, 1)
    A = 0.0; Sy = 0.0; Sz = 0.0
    for i in 1:nelem
        n1 = Int(cs_elem[i, 2]); n2 = Int(cs_elem[i, 3]); t = cs_elem[i, 4]
        y1 = cs_node[n1, 2]; z1 = cs_node[n1, 3]
        y2 = cs_node[n2, 2]; z2 = cs_node[n2, 3]
        L = sqrt((y2 - y1)^2 + (z2 - z1)^2)
        dA = t * L
        A  += dA
        Sy += dA * (z1 + z2) / 2.0
        Sz += dA * (y1 + y2) / 2.0
    end
    yCG = Sz / A; zCG = Sy / A

    Iyr = 0.0; Izr = 0.0; Iyzr = 0.0
    for i in 1:nelem
        n1 = Int(cs_elem[i, 2]); n2 = Int(cs_elem[i, 3]); t = cs_elem[i, 4]
        y1 = cs_node[n1, 2]; z1 = cs_node[n1, 3]
        y2 = cs_node[n2, 2]; z2 = cs_node[n2, 3]
        L = sqrt((y2 - y1)^2 + (z2 - z1)^2)
        dIyr  = t * L * (z1^2 + z1*z2 + z2^2) / 3.0
        dIzr  = t * L * (y1^2 + y1*y2 + y2^2) / 3.0
        dIyzr = t * L * (2.0*y1*z1 + y1*z2 + y2*z1 + 2.0*y2*z2) / 6.0
        Iyr  += dIyr; Izr  += dIzr; Iyzr += dIyzr
    end
    Iy = Iyr - A * zCG^2
    Iz = Izr - A * yCG^2
    return A, yCG, zCG, Iyr, Iy, Izr, Iz
end

function cs_data(h, b, c, tb, th, tc)
    cs_node, cs_elem = cs_def(h, b, c, tb, th, tc)
    return cs_prop(cs_node, cs_elem)
end

# ─────────────────────────────────────────────────────────────────────────────
# GAUSS INTEGRATION
# ─────────────────────────────────────────────────────────────────────────────

function GL3(np::Int)
    if np == 1
        locx = [1.0/3.0]; locy = [1.0/3.0]; wei = [1.0]
    elseif np == 3
        locx = [1.0/6.0, 2.0/3.0, 1.0/6.0]
        locy = [1.0/6.0, 1.0/6.0, 2.0/3.0]
        wei  = [1.0/3.0, 1.0/3.0, 1.0/3.0]
    elseif np == 4
        locx = [1.0/3.0, 0.6, 0.2, 0.2]
        locy = [1.0/3.0, 0.2, 0.6, 0.2]
        wei  = [-27.0/48.0, 25.0/48.0, 25.0/48.0, 25.0/48.0]
    else
        error("GL3: np must be 1, 3, or 4")
    end
    return locx, locy, wei
end

function GL4(np::Int)
    if np == 1
        locx = [0.0]; locy = [0.0]; wei = [4.0]
    elseif np == 4
        a = 1.0 / sqrt(3.0)
        locx = [-a,  a,  a, -a]
        locy = [-a, -a,  a,  a]
        wei  = [1.0, 1.0, 1.0, 1.0]
    elseif np == 9
        a = sqrt(0.6)
        locx = [-a, 0.0, a, -a, 0.0, a, -a, 0.0, a]
        locy = [-a, -a, -a, 0.0, 0.0, 0.0, a, a, a]
        w1 = 5.0/9.0; w2 = 8.0/9.0; w3 = 5.0/9.0
        wvec = [w1, w2, w3]
        wei = zeros(9)
        idx = 1
        for j in 1:3, i in 1:3
            wei[idx] = wvec[i] * wvec[j]
            idx += 1
        end
    else
        error("GL4: np must be 1, 4, or 9")
    end
    return locx, locy, wei
end

# ─────────────────────────────────────────────────────────────────────────────
# FINITE ELEMENT DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────

function FE_props()
    FE_dat = Dict{Int, FEData}()

    # t3
    FE_dat[11] = FEData(3, 3, 6, 6, false, false, 0, "Mindlin", 9, 9, false, false, 0, 100, 18, 1)

    # t31
    FE_dat[12] = FEData(3, 3, 8, 6, true, false, 0, "Mindlin", 12, 9, true, false, 0, 100, 18, 1)

    # t33
    FE_dat[13] = FEData(3, 3, 12, 6, true, false, 0, "Mindlin", 18, 9, true, false, 0, 100, 18, 1)

    # t6
    FE_dat[14] = FEData(3, 6, 12, 12, false, true, 4, "Mindlin", 18, 18, false, true, 4, 100, 18, 4)

    # q4
    FE_dat[21] = FEData(4, 4, 8, 8, false, true, 4, "Mindlin", 12, 12, false, true, 4, 100, 24, 4)

    # q42
    FE_dat[22] = FEData(4, 4, 12, 8, true, true, 4, "Mindlin", 18, 12, true, true, 9, 100, 24, 4)

    # q8
    FE_dat[23] = FEData(4, 8, 16, 16, false, true, 9, "Mindlin", 24, 24, false, true, 9, 100, 24, 9)

    # q81
    FE_dat[24] = FEData(4, 8, 18, 16, true, true, 9, "Mindlin", 27, 24, true, true, 9, 100, 24, 9)

    return FE_dat
end

# ─────────────────────────────────────────────────────────────────────────────
# COORDINATE TRANSFORMATIONS
# ─────────────────────────────────────────────────────────────────────────────

function rotate3d(T3::AbstractMatrix, nodenr::Int, dpn::Int)
    tot = nodenr * dpn
    T = zeros(tot, tot)
    for inod in 1:nodenr
        i0 = (inod - 1) * dpn
        T[i0+1:i0+3, i0+1:i0+3] .= T3
        if dpn == 6
            T[i0+4:i0+6, i0+4:i0+6] .= T3
        end
    end
    return T
end

function ct_3node_g2e(elems_row::AbstractVector, nodes::AbstractMatrix)
    n1 = Int(elems_row[1]); n2 = Int(elems_row[2]); n3 = Int(elems_row[3])
    P1g = nodes[n1, 1:3]; P2g = nodes[n2, 1:3]; P3g = nodes[n3, 1:3]
    v21 = P2g .- P1g
    j1 = v21 ./ norm(v21)
    v31 = P3g .- P1g
    c = cross(v21, v31)
    j3 = c ./ norm(c)
    j2 = cross(j3, j1)
    T = hcat(j1, j2, j3)
    P1e = T' * (P1g .- P1g)
    P2e = T' * (P2g .- P1g)
    P3e = T' * (P3g .- P1g)
    return T, P1e, P2e, P3e
end

function ct_4node_g2e(elems_row::AbstractVector, nodes::AbstractMatrix)
    n1 = Int(elems_row[1]); n2 = Int(elems_row[2])
    n3 = Int(elems_row[3]); n4 = Int(elems_row[4])
    P1g = nodes[n1, 1:3]; P2g = nodes[n2, 1:3]
    P3g = nodes[n3, 1:3]; P4g = nodes[n4, 1:3]
    v21 = P2g .- P1g
    v41 = P4g .- P1g
    c = cross(v21, v41)
    j3 = c ./ norm(c)
    P13 = (P1g .+ P3g) ./ 2.0
    P24 = (P2g .+ P4g) ./ 2.0
    v13 = P3g .- P1g; v24 = P4g .- P2g
    v13 = v13 ./ norm(v13); v24 = v24 ./ norm(v24)
    j1 = (v13 .- v24) ./ norm(v13 .- v24)
    j2 = (v13 .+ v24) ./ norm(v13 .+ v24)
    T = hcat(j1, j2, j3)
    P1e = T' * (P1g .- P13)
    P2e = T' * (P2g .- P13)
    P3e = T' * (P3g .- P13)
    P4e = T' * (P4g .- P13)
    return T, P1e, P2e, P3e, P4e
end

function plate_transform(P1::AbstractVector, P2::AbstractVector,
                         P3::AbstractVector, P4::AbstractVector)
    v21 = P2 .- P1; v41 = P4 .- P1
    c = cross(v21, v41)
    j3 = c ./ norm(c)
    j1 = v21 ./ norm(v21)
    j2 = cross(j3, j1)
    T = hcat(j1, j2, j3)
    P0 = copy(P1)
    return P0, T
end

function dir_convert(dir_char)
    d = string(dir_char)
    if uppercase(d) in ("X", "RX")
        return 1
    elseif uppercase(d) in ("Y", "RY")
        return 2
    elseif uppercase(d) in ("Z", "RZ")
        return 3
    else
        return 1
    end
end

function collect_pl_data(PLATE_table)
    npl = length(PLATE_table)
    Pplg = zeros(npl, 5, 3)
    Ppls = zeros(npl, 5, 3)
    Tpl  = zeros(npl, 3, 3)
    for ipl in 1:npl
        Pg1 = PLATE_table[ipl][1:3]
        Pg2 = PLATE_table[ipl][4:6]
        Pg3 = PLATE_table[ipl][7:9]
        Pg4 = PLATE_table[ipl][10:12]
        P0, T = plate_transform(Pg1, Pg2, Pg3, Pg4)
        Ps1 = T' * (Pg1 .- Pg1)
        Ps2 = T' * (Pg2 .- Pg1)
        Ps3 = T' * (Pg3 .- Pg1)
        Ps4 = T' * (Pg4 .- Pg1)
        Ps0 = T' * (P0 .- Pg1)
        Pplg[ipl, 1, :] = Pg1; Pplg[ipl, 2, :] = Pg2; Pplg[ipl, 3, :] = Pg3; Pplg[ipl, 4, :] = Pg4; Pplg[ipl, 5, :] = P0
        Ppls[ipl, 1, :] = Ps1; Ppls[ipl, 2, :] = Ps2; Ppls[ipl, 3, :] = Ps3; Ppls[ipl, 4, :] = Ps4; Ppls[ipl, 5, :] = Ps0
        Tpl[ipl, :, :] = T
    end
    return Pplg, Ppls, Tpl
end

# ─────────────────────────────────────────────────────────────────────────────
# MESH GENERATION (RECTANGULAR & DISTMESH)
# ─────────────────────────────────────────────────────────────────────────────

function collect_anchor(PLATE_table, SUPPORT_table, SPR_EXT_table, tol)
    nplate = length(PLATE_table)
    p_anchor = [Matrix{Float64}(undef, 0, 3) for _ in 1:nplate]

    for ipl in 1:nplate
        P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
        P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
        P0, T = plate_transform(P1, P2, P3, P4)
        Q = zeros(4, 3)
        Q[1, :] = (P1 .- P0)' * T
        Q[2, :] = (P2 .- P0)' * T
        Q[3, :] = (P3 .- P0)' * T
        Q[4, :] = (P4 .- P0)' * T

        for iside in 1:4
            i1 = iside; i2 = (iside == 4) ? 1 : iside + 1
            Q1 = Q[i1, :]; Q2 = Q[i2, :]

            # Loop over SPR_EXT_table (or CONDSUP_table)
            for isup in 1:length(SPR_EXT_table)
                entry = SPR_EXT_table[isup]
                if isempty(entry); continue; end
                stype = entry[1]
                if stype == "point"
                    jpl = entry[2]
                    Pp = zeros(3)
                    Pp[1] = Float64(entry[3]); Pp[2] = Float64(entry[4]); Pp[3] = 0.0
                    P1_j = PLATE_table[jpl][1:3]; P2_j = PLATE_table[jpl][4:6]
                    P3_j = PLATE_table[jpl][7:9]; P4_j = PLATE_table[jpl][10:12]
                    P0_j, T_j = plate_transform(P1_j, P2_j, P3_j, P4_j)
                    Pg = T_j * Pp .+ P0_j
                    Ps = T' * (Pg .- P0)
                    if abs((Ps[1]-Q2[1])*(Q1[2]-Q2[2]) - (Q1[1]-Q2[1])*(Ps[2]-Q2[2])) < tol &&
                       Ps[1] >= min(Q1[1],Q2[1])-tol && Ps[1] <= max(Q1[1],Q2[1])+tol &&
                       Ps[2] >= min(Q1[2],Q2[2])-tol && Ps[2] <= max(Q1[2],Q2[2])+tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], reshape(Pg, 1, 3))
                    end
                elseif stype == "surf"
                    jpl = entry[2]
                    Pp = zeros(3)
                    Pp[1] = Float64(entry[3]); Pp[2] = Float64(entry[4]); Pp[3] = 0.0
                    P1_j = PLATE_table[jpl][1:3]; P2_j = PLATE_table[jpl][4:6]
                    P3_j = PLATE_table[jpl][7:9]; P4_j = PLATE_table[jpl][10:12]
                    P0_j, T_j = plate_transform(P1_j, P2_j, P3_j, P4_j)
                    Pg = T_j * Pp .+ P0_j
                    Ps = T' * (Pg .- P0)
                    if abs((Ps[1]-Q2[1])*(Q1[2]-Q2[2]) - (Q1[1]-Q2[1])*(Ps[2]-Q2[2])) < tol &&
                       Ps[1] >= min(Q1[1],Q2[1])-tol && Ps[1] <= max(Q1[1],Q2[1])+tol &&
                       Ps[2] >= min(Q1[2],Q2[2])-tol && Ps[2] <= max(Q1[2],Q2[2])+tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], reshape(Pg, 1, 3))
                    end
                end
            end

            # Loop over SUPPORT_table
            for isup in 1:length(SUPPORT_table)
                entry = SUPPORT_table[isup]
                if isempty(entry); continue; end
                stype = entry[1]
                if stype == "point"
                    jpl = entry[2]
                    Pp = zeros(3)
                    Pp[1] = Float64(entry[3]); Pp[2] = Float64(entry[4]); Pp[3] = 0.0
                    P1_j = PLATE_table[jpl][1:3]; P2_j = PLATE_table[jpl][4:6]
                    P3_j = PLATE_table[jpl][7:9]; P4_j = PLATE_table[jpl][10:12]
                    P0_j, T_j = plate_transform(P1_j, P2_j, P3_j, P4_j)
                    Pg = T_j * Pp .+ P0_j
                    Ps = T' * (Pg .- P0)
                    if abs((Ps[1]-Q2[1])*(Q1[2]-Q2[2]) - (Q1[1]-Q2[1])*(Ps[2]-Q2[2])) < tol &&
                       Ps[1] >= min(Q1[1],Q2[1])-tol && Ps[1] <= max(Q1[1],Q2[1])+tol &&
                       Ps[2] >= min(Q1[2],Q2[2])-tol && Ps[2] <= max(Q1[2],Q2[2])+tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], reshape(Pg, 1, 3))
                    end
                end
            end
        end
    end
    return p_anchor
end

function mesh4_rect(P1, P2, P3, P4, p_anchor_mat, size_x, size_y, tol)
    P0, T = plate_transform(P1, P2, P3, P4)
    P1s = T' * (P1 .- P0)
    P3s = T' * (P3 .- P0)

    x_anch = [P1s[1], P3s[1]]
    y_anch = [P1s[2], P3s[2]]

    for i in 1:size(p_anchor_mat, 1)
        Ps = T' * (p_anchor_mat[i, :] .- P0)
        if minimum(abs.(x_anch .- Ps[1])) > tol
            push!(x_anch, Ps[1])
        end
        if minimum(abs.(y_anch .- Ps[2])) > tol
            push!(y_anch, Ps[2])
        end
    end

    sort!(x_anch); sort!(y_anch)
    x_nod = [x_anch[1]]; y_nod = [y_anch[1]]

    for i in 1:(length(x_anch) - 1)
        L = x_anch[i+1] - x_anch[i]
        n = max(1, round(Int, L / size_x))
        dL = L / n
        for j in 1:n
            push!(x_nod, x_anch[i] + j * dL)
        end
    end

    for i in 1:(length(y_anch) - 1)
        L = y_anch[i+1] - y_anch[i]
        n = max(1, round(Int, L / size_y))
        dL = L / n
        for j in 1:n
            push!(y_nod, y_anch[i] + j * dL)
        end
    end

    nx = length(x_nod); ny = length(y_nod)
    nodes_loc = zeros(nx * ny, 3)
    idx = 1
    for iy in 1:ny, ix in 1:nx
        nodes_loc[idx, 1] = x_nod[ix]
        nodes_loc[idx, 2] = y_nod[iy]
        idx += 1
    end

    elems = zeros(Int, (nx - 1) * (ny - 1), 4)
    ie = 1
    for iy in 1:(ny - 1), ix in 1:(nx - 1)
        n1 = (iy - 1) * nx + ix
        n2 = n1 + 1
        n3 = iy * nx + ix + 1
        n4 = iy * nx + ix
        elems[ie, :] = [n1, n2, n3, n4]
        ie += 1
    end

    nodes = zeros(nx * ny, 3)
    for inod in 1:(nx * ny)
        nodes[inod, :] = T * nodes_loc[inod, :] .+ P0
    end

    return nodes, elems
end

function edges_by_plate(pl_elems)
    npl = length(pl_elems)
    pl_el_edges = Vector{Matrix{Int}}(undef, npl)
    pl_edge_nodes = Vector{Matrix{Int}}(undef, npl)

    for ipl in 1:npl
        ele = size(pl_elems[ipl], 1)
        if ele == 0
            pl_el_edges[ipl] = zeros(Int, 0, 4)
            pl_edge_nodes[ipl] = zeros(Int, 0, 4)
            continue
        end

        el_edge = zeros(Int, ele, 4)
        edge_nodes = zeros(Int, ele * 4, 2)
        iedge = 0

        for iel in 1:ele
            el_nod = pl_elems[ipl][iel, :]
            n_node = (el_nod[4] == 0) ? 3 : 4
            for inod in 1:n_node
                inod1 = inod
                inod2 = (inod == n_node) ? 1 : inod + 1
                nod1 = el_nod[inod1]; nod2 = el_nod[inod2]
                n1 = min(nod1, nod2); n2 = max(nod1, nod2)

                matched = false
                for k in 1:iedge
                    if edge_nodes[k, 1] == n1 && edge_nodes[k, 2] == n2
                        el_edge[iel, inod] = k
                        matched = true
                        break
                    end
                end
                if !matched
                    iedge += 1
                    edge_nodes[iedge, 1] = n1
                    edge_nodes[iedge, 2] = n2
                    el_edge[iel, inod] = iedge
                end
            end
        end

        pl_el_edges[ipl] = el_edge
        edge_nodes_out = zeros(Int, iedge, 4)
        edge_nodes_out[:, 1:2] .= edge_nodes[1:iedge, :]
        pl_edge_nodes[ipl] = edge_nodes_out
    end

    return pl_el_edges, pl_edge_nodes
end

# ─────────────────────────────────────────────────────────────────────────────
# 2D DELAUNAY TRIANGULATION & DISTMESH
# ─────────────────────────────────────────────────────────────────────────────

function delaunay_2d(points::Matrix{Float64})
    n = size(points, 1)
    if n < 3
        return zeros(Int, 0, 3)
    end

    xmin = minimum(points[:, 1]); xmax = maximum(points[:, 1])
    ymin = minimum(points[:, 2]); ymax = maximum(points[:, 2])
    dx = xmax - xmin; dy = ymax - ymin
    dmax = max(dx, dy, 1e-6)
    xmid = (xmin + xmax) / 2.0; ymid = (ymin + ymax) / 2.0

    p_super1 = [xmid - 20.0 * dmax, ymid - dmax]
    p_super2 = [xmid, ymid + 20.0 * dmax]
    p_super3 = [xmid + 20.0 * dmax, ymid - dmax]

    pts = Matrix{Float64}(undef, n + 3, 2)
    pts[1:n, :] = points
    pts[n+1, :] = p_super1
    pts[n+2, :] = p_super2
    pts[n+3, :] = p_super3

    triangles = Set{NTuple{3, Int}}()
    push!(triangles, (n + 1, n + 2, n + 3))

    function in_circle(i::Int, j::Int, k::Int, px::Float64, py::Float64)
        ax = pts[i, 1] - px; ay = pts[i, 2] - py
        bx = pts[j, 1] - px; by = pts[j, 2] - py
        cx = pts[k, 1] - px; cy = pts[k, 2] - py

        det = (ax^2 + ay^2) * (bx * cy - cx * by) -
              (bx^2 + by^2) * (ax * cy - cx * ay) +
              (cx^2 + cy^2) * (ax * by - bx * ay)

        orient = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
        return orient > 0 ? det > 1e-14 : det < -1e-14
    end

    for p_idx in 1:n
        px = pts[p_idx, 1]; py = pts[p_idx, 2]

        bad_triangles = NTuple{3, Int}[]
        for tri in triangles
            if in_circle(tri[1], tri[2], tri[3], px, py)
                push!(bad_triangles, tri)
            end
        end

        edge_counts = Dict{Tuple{Int, Int}, Int}()
        for tri in bad_triangles
            edges = [
                (min(tri[1], tri[2]), max(tri[1], tri[2])),
                (min(tri[2], tri[3]), max(tri[2], tri[3])),
                (min(tri[3], tri[1]), max(tri[3], tri[1]))
            ]
            for e in edges
                edge_counts[e] = get(edge_counts, e, 0) + 1
            end
            delete!(triangles, tri)
        end

        for (e, count) in edge_counts
            if count == 1
                ax = pts[e[1], 1]; ay = pts[e[1], 2]
                bx = pts[e[2], 1]; by = pts[e[2], 2]
                orient = (bx - ax) * (py - ay) - (by - ay) * (px - ax)
                if orient > 0
                    push!(triangles, (e[1], e[2], p_idx))
                else
                    push!(triangles, (e[2], e[1], p_idx))
                end
            end
        end
    end

    valid_tris = NTuple{3, Int}[]
    for tri in triangles
        if tri[1] <= n && tri[2] <= n && tri[3] <= n
            push!(valid_tris, tri)
        end
    end

    T = Matrix{Int}(undef, length(valid_tris), 3)
    for (idx, tri) in enumerate(valid_tris)
        T[idx, 1] = tri[1]
        T[idx, 2] = tri[2]
        T[idx, 3] = tri[3]
    end
    return T
end

function deduplicate_points(p::Matrix{Float64}, tol::Float64=1e-6)
    n = size(p, 1)
    keep = trues(n)
    for i in 1:n
        if !keep[i]; continue; end
        for j in i+1:n
            if !keep[j]; continue; end
            if abs(p[i, 1] - p[j, 1]) < tol && abs(p[i, 2] - p[j, 2]) < tol
                keep[j] = false
            end
        end
    end
    return p[keep, :]
end

function distmesh(fd, fh, h0::Float64, bbox::Matrix{Float64}, p_fix::Matrix{Float64}=zeros(0, 2); it_max::Int=200)
    dp_tol = -0.001 * h0
    dtrm_tol = -0.001 * h0
    rt_tol = 0.3
    F_scale = 1.2
    dp_scale = 0.2
    dpc_tol = 0.001 * h0
    gradeps = sqrt(eps(Float64)) * h0

    x_range = bbox[1, 1]:h0:bbox[2, 1]
    y_range = bbox[1, 2]:(h0 * sqrt(3.0)/2.0):bbox[2, 2]

    pts_list = Vector{Float64}[]
    for (iy, y) in enumerate(y_range)
        shift = (iy % 2 == 0) ? h0 / 2.0 : 0.0
        for x in x_range
            push!(pts_list, [x + shift, y])
        end
    end
    p = reduce(vcat, [reshape(pt, 1, 2) for pt in pts_list])

    d = fd(p)
    p = p[d .< -dp_tol, :]

    if size(p_fix, 1) > 0
        p_fix = deduplicate_points(p_fix, 1e-6)
        n_fix = size(p_fix, 1)
        keep = trues(size(p, 1))
        for i in 1:size(p, 1)
            for j in 1:n_fix
                if norm(p[i, :] - p_fix[j, :]) < 0.5 * h0
                    keep[i] = false
                    break
                end
            end
        end
        p = [p_fix; p[keep, :]]
    else
        n_fix = 0
    end
    p = deduplicate_points(p, 1e-6)
    n_p = size(p, 1)

    p0 = fill(Inf, n_p, 2)
    t = zeros(Int, 0, 3)

    for it in 1:it_max
        delta_p_max = maximum(sqrt.(sum((p .- p0).^2, dims=2)))
        if rt_tol * h0 < delta_p_max
            p0 = copy(p)
            t = delaunay_2d(p)
            pc = (p[t[:, 1], :] + p[t[:, 2], :] + p[t[:, 3], :]) ./ 3.0
            dist_c = fd(pc)
            valid = dist_c .< dtrm_tol
            t = t[valid, :]

            for i in 1:size(t, 1)
                ax = p[t[i, 2], 1] - p[t[i, 1], 1]; ay = p[t[i, 2], 2] - p[t[i, 1], 2]
                bx = p[t[i, 3], 1] - p[t[i, 1], 1]; by = p[t[i, 3], 2] - p[t[i, 1], 2]
                if ax * by - ay * bx < 0
                    t[i, 1], t[i, 2] = t[i, 2], t[i, 1]
                end
            end
        end

        if size(t, 1) == 0; break; end

        edge_set = Set{Tuple{Int, Int}}()
        for i in 1:size(t, 1)
            push!(edge_set, (min(t[i, 1], t[i, 2]), max(t[i, 1], t[i, 2])))
            push!(edge_set, (min(t[i, 2], t[i, 3]), max(t[i, 2], t[i, 3])))
            push!(edge_set, (min(t[i, 3], t[i, 1]), max(t[i, 3], t[i, 1])))
        end
        edges = collect(edge_set)
        n_edges = length(edges)

        bars = zeros(n_edges, 2)
        L = zeros(n_edges)
        for (ie, e) in enumerate(edges)
            diff = p[e[1], :] - p[e[2], :]
            bars[ie, :] = diff
            L[ie] = norm(diff)
        end

        midpts = zeros(n_edges, 2)
        for (ie, e) in enumerate(edges)
            midpts[ie, :] = 0.5 * (p[e[1], :] + p[e[2], :])
        end
        hbars = fh(midpts)
        L_target = hbars .* F_scale * (sum(L.^2) / sum(hbars.^2))^0.5

        F = max.(L_target .- L, 0.0)
        F_bar = (F ./ max.(L, 1e-12)) .* bars

        delta_p = zeros(size(p, 1), 2)
        for (ie, e) in enumerate(edges)
            delta_p[e[1], :] .+= F_bar[ie, :]
            delta_p[e[2], :] .-= F_bar[ie, :]
        end

        if n_fix > 0
            delta_p[1:n_fix, :] .= 0.0
        end
        delta_p .*= dp_scale
        p .+= delta_p

        d = fd(p)
        ix = findall(d .> 0)
        filter!(idx -> idx > n_fix, ix)
        if !isempty(ix)
            for idx in ix
                p_cur = p[idx, :]
                d_val = d[idx]
                d_x = fd(reshape(p_cur + [gradeps, 0.0], 1, 2))[1]
                d_y = fd(reshape(p_cur + [0.0, gradeps], 1, 2))[1]
                grad = [(d_x - d_val) / gradeps, (d_y - d_val) / gradeps]
                gnm = norm(grad)^2
                if gnm > 1e-12
                    p[idx, :] .-= (d_val / gnm) .* grad
                end
            end
        end

        max_dp = maximum(sqrt.(sum(delta_p.^2, dims=2)))
        if max_dp < dpc_tol
            break
        end
    end

    t = delaunay_2d(p)
    pc = (p[t[:, 1], :] + p[t[:, 2], :] + p[t[:, 3], :]) ./ 3.0
    dist_c = fd(pc)
    t = t[dist_c .< dtrm_tol, :]
    for i in 1:size(t, 1)
        ax = p[t[i, 2], 1] - p[t[i, 1], 1]; ay = p[t[i, 2], 2] - p[t[i, 1], 2]
        bx = p[t[i, 3], 1] - p[t[i, 1], 1]; by = p[t[i, 3], 2] - p[t[i, 1], 2]
        if ax * by - ay * bx < 0
            t[i, 1], t[i, 2] = t[i, 2], t[i, 1]
        end
    end

    return p, t
end

# ─────────────────────────────────────────────────────────────────────────────
# MODEL GENERATION (DISTMESH VERSION - Sandor Adany, latest)
# ─────────────────────────────────────────────────────────────────────────────

function main_modelgen_dm(PLATE_table, LOAD_table, SUPPORT_table, SPR_EXT_table, HOLE_table,
                          MPC_SECTION_table, MPC_HOLE_table, FE_dat, fe_size, mesh_typ, rat_q2t, tolers)
    tol = tolers.zero_dist
    p_anchor = collect_anchor(PLATE_table, SUPPORT_table, SPR_EXT_table, tol)

    pl_nodes = Vector{Matrix{Float64}}(undef, length(PLATE_table))
    pl_elems = Vector{Matrix{Int}}(undef, length(PLATE_table))

    coord0 = zeros(0, 2)

    # Mesh plates 2:nplate (the upright column walls)
    for ipl in 2:length(PLATE_table)
        P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
        P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
        nod, ele = mesh4_rect(P1, P2, P3, P4, p_anchor[ipl], fe_size[1], fe_size[2], tol)
        pl_nodes[ipl] = nod
        nel = size(ele, 1)
        pl_elems_ipl = zeros(Int, nel, 9)
        pl_elems_ipl[:, 1:4] = ele
        pl_elems[ipl] = pl_elems_ipl

        nod0_idx = findall(abs.(nod[:, 3]) .< tol)
        if !isempty(nod0_idx)
            coord0 = vcat(coord0, nod[nod0_idx, 1:2])
        end
    end

    # Mesh plate 1 (the base plate) using distmesh
    ipl = 1
    P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
    P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]

    p_up = coord0
    p_hole_mid = zeros(0, 2)
    p_hole = zeros(0, 2)
    if length(HOLE_table) >= 1 && HOLE_table[1] !== nothing
        for c in HOLE_table[1].coord
            p_hole_mid = vcat(p_hole_mid, reshape(c[1:2], 1, 2))
        end
        for h in HOLE_table[1].holes
            p_hole = vcat(p_hole, h[:, 1:2])
        end
    end
    p_corner = [P1[1:2]'; P2[1:2]'; P3[1:2]'; P4[1:2]']
    pp = vcat(p_corner, p_up, p_hole, p_hole_mid)

    a0 = fe_size[1]
    fd = p -> -min.(min.(min.(-P1[2] .+ p[:, 2], P3[2] .- p[:, 2]), -P1[1] .+ p[:, 1]), P2[1] .- p[:, 1])
    fh = p -> ones(size(p, 1))

    p, t = distmesh(fd, fh, a0, [P1[1] P1[2]; P3[1] P3[2]], pp)

    pl_elems_1 = zeros(Int, size(t, 1), 9)
    pl_elems_1[:, 1:3] = t
    pl_elems[1] = pl_elems_1
    pl_nodes[1] = hcat(p, zeros(size(p, 1), 1))

    pl_el_edges, pl_edge_nodes = edges_by_plate(pl_elems)

    # Merge elements and nodes
    tot_ele = sum(size(pl_elems[i], 1) for i in 1:length(PLATE_table))
    tot_nod = sum(size(pl_nodes[i], 1) for i in 1:length(PLATE_table))

    el_nodes = zeros(Int, tot_ele, 9)
    el_props = zeros(tot_ele, 5)
    nodes = zeros(tot_nod, 5)

    ele0 = 0; nod0 = 0
    for ipl in 1:length(PLATE_table)
        ele = size(pl_elems[ipl], 1)
        nod = size(pl_nodes[ipl], 1)
        for iel in 1:ele
            row = pl_elems[ipl][iel, :]
            nz = findall(row .> 0)
            el_nodes[ele0 + iel, nz] = row[nz] .+ nod0
        end
        el_props[ele0+1:ele0+ele, 1] .= ipl
        pl_edge_nodes[ipl][:, 3:4] .+= nod0

        # Element FE type assignment
        ind_tri = findall(pl_elems[ipl][:, 4] .== 0)
        el_props[ind_tri .+ ele0, 2] .= mesh_typ[ipl, 1]
        ind_quad = findall(pl_elems[ipl][:, 4] .> 0)
        el_props[ind_quad .+ ele0, 2] .= mesh_typ[ipl, 2]

        el_props[ele0+1:ele0+ele, 3] .= PLATE_table[ipl][13]
        el_props[ele0+1:ele0+ele, 4] .= PLATE_table[ipl][14]

        nodes[nod0+1:nod0+nod, 1:3] = pl_nodes[ipl]
        nodes[nod0+1:nod0+nod, 4] = (nod0+1):(nod0+nod)

        ele0 += ele; nod0 += nod
    end

    # Add hole information
    el_props[:, 5] .= 1.0  # default: no hole
    for ipl in 1:length(PLATE_table)
        if ipl <= length(HOLE_table) && HOLE_table[ipl] !== nothing
            C1 = PLATE_table[ipl][1:3]; C2 = PLATE_table[ipl][4:6]
            C3 = PLATE_table[ipl][7:9]; C4 = PLATE_table[ipl][10:12]
            C0, T = plate_transform(C1, C2, C3, C4)
            for ihole in 1:length(HOLE_table[ipl].holes)
                x0 = HOLE_table[ipl].coord[ihole][1]
                y0 = HOLE_table[ipl].coord[ihole][2]
                R  = HOLE_table[ipl].radius[ihole]
                sec = HOLE_table[ipl].sector[ihole]
                R_ = R * cos(pi / sec)
                C = T * [x0, y0, 0.0] .+ C1
                for iel in 1:size(el_nodes, 1)
                    fe_id = Int(el_props[iel, 2])
                    elemnodenr = FE_dat[fe_id].nodenr
                    midpoint = vec(mean(nodes[el_nodes[iel, 1:elemnodenr], 1:3], dims=1))
                    if norm(midpoint .- C) < R_ + tol
                        el_props[iel, 5] = -1.0
                        el_props[iel, 3] /= 100.0
                    end
                end
            end
        end
    end

    # Deduplicate nodes across plate interfaces
    for in1 in 1:(size(nodes, 1) - 1)
        if nodes[in1, 4] == in1
            p1 = nodes[in1, 1:3]
            for in2 in (in1 + 1):size(nodes, 1)
                if abs(nodes[in2, 1] - p1[1]) < tol &&
                   abs(nodes[in2, 2] - p1[2]) < tol &&
                   abs(nodes[in2, 3] - p1[3]) < tol
                    nodes[in2, 4] = 0.0
                    nodes[in2, 5] = in1
                end
            end
        end
    end

    # Renumber active nodes
    nr = 0
    for i in 1:size(nodes, 1)
        if nodes[i, 4] > 0
            nr += 1
            nodes[i, 5] = nr
        end
    end
    for i in 1:size(nodes, 1)
        if nodes[i, 4] == 0
            tgt = Int(nodes[i, 5])
            nodes[i, 5] = nodes[tgt, 5]
        end
    end
    nodenr = nr

    # Add MPC beam element data
    n_mpc_sec = (size(MPC_SECTION_table, 1) > 0) ? size(MPC_SECTION_table, 1) : 0
    n_mpc_hol = (size(MPC_HOLE_table, 1) > 0) ? size(MPC_HOLE_table, 1) : 0

    beam_nodes = zeros(Int, 0, 3)
    beam_props = zeros(Float64, 0, 7)

    # MPC for cross-sections
    for i in 1:n_mpc_sec
        nodnr = size(nodes, 1) + 1
        H_sec = MPC_SECTION_table[i, 1]
        nods = findall((abs.(nodes[:, 3] .- H_sec) .< tol) .&& (nodes[:, 4] .> 0))
        for inod in nods
            beam_nodes = vcat(beam_nodes, [nodnr inod 0])
            beam_props = vcat(beam_props, reshape([1.0, 1.0, MPC_SECTION_table[i, 5:9]...], 1, 7))
        end
        nodenr += 1
        new_row = [MPC_SECTION_table[i, 2:4]... nodnr nodenr]
        nodes = vcat(nodes, reshape(new_row, 1, 5))
    end

    # MPC for bolt holes
    for i in 1:n_mpc_hol
        ipl_h = Int(MPC_HOLE_table[i, 1])
        iho = Int(MPC_HOLE_table[i, 2])
        C1 = PLATE_table[ipl_h][1:3]; C2 = PLATE_table[ipl_h][4:6]
        C3 = PLATE_table[ipl_h][7:9]; C4 = PLATE_table[ipl_h][10:12]
        C0, T = plate_transform(C1, C2, C3, C4)
        sec = HOLE_table[ipl_h].sector[iho]
        x0 = HOLE_table[ipl_h].coord[iho][1]
        y0 = HOLE_table[ipl_h].coord[iho][2]
        C = T * [x0, y0, 0.0] .+ C1

        nod_m_list = findall((abs.(nodes[:, 1] .- C[1]) .< tol) .&&
                             (abs.(nodes[:, 2] .- C[2]) .< tol) .&&
                             (abs.(nodes[:, 3] .- C[3]) .< tol) .&&
                             (nodes[:, 4] .> 0))
        if !isempty(nod_m_list)
            nod_m = nod_m_list[1]
            for inod in 1:sec
                coord = [HOLE_table[ipl_h].holes[iho][inod, 1:2]..., 0.0]
                C_p = T * coord .+ C1
                nod_p_list = findall((abs.(nodes[:, 1] .- C_p[1]) .< tol) .&&
                                     (abs.(nodes[:, 2] .- C_p[2]) .< tol) .&&
                                     (abs.(nodes[:, 3] .- C_p[3]) .< tol) .&&
                                     (nodes[:, 4] .> 0))
                if !isempty(nod_p_list)
                    nod_p = nod_p_list[1]
                    beam_nodes = vcat(beam_nodes, [nod_m nod_p 0])
                    beam_props = vcat(beam_props, reshape([1.0, 1.0, MPC_HOLE_table[i, 3:7]...], 1, 7))
                end
            end
        end
    end

    return nodes, el_nodes, el_props, nodenr, pl_el_edges, pl_edge_nodes, beam_nodes, beam_props
end

# ─────────────────────────────────────────────────────────────────────────────
# BEAM FINITE ELEMENTS (EULER-BERNOULLI 3D BEAMS)
# ─────────────────────────────────────────────────────────────────────────────

function ke_beam_EB(P1g::AbstractVector, P2g::AbstractVector, P3g::AbstractVector,
                    A::Real, Iy::Real, Iz::Real, J::Real, E::Real, G::Real)
    L = norm(P2g .- P1g)
    norm_vec = cross(P2g .- P1g, P3g .- P1g)
    j3 = norm_vec ./ norm(norm_vec)
    j1 = (P2g .- P1g) ./ L
    j2 = cross(j3, j1)
    T = hcat(j1, j2, j3)

    kx   = E * A / L
    kxx  = G * J / L^3
    ky2  = 2.0 * E * Iy / L; ky4  = 4.0 * E * Iy / L; ky6  = 6.0 * E * Iy / L^2; ky12 = 12.0 * E * Iy / L^3
    kz2  = 2.0 * E * Iz / L; kz4  = 4.0 * E * Iz / L; kz6  = 6.0 * E * Iz / L^2; kz12 = 12.0 * E * Iz / L^3

    ke = [
        kx     0.0   0.0   0.0   0.0   0.0  -kx     0.0   0.0   0.0   0.0   0.0
        0.0   kz12   0.0   0.0   0.0   kz6   0.0  -kz12   0.0   0.0   0.0   kz6
        0.0    0.0  ky12   0.0  -ky6   0.0   0.0    0.0 -ky12   0.0  -ky6   0.0
        0.0    0.0   0.0   kxx   0.0   0.0   0.0    0.0   0.0  -kxx   0.0   0.0
        0.0    0.0  -ky6   0.0   ky4   0.0   0.0    0.0   ky6   0.0   ky2   0.0
        0.0   kz6    0.0   0.0   0.0   kz4   0.0  -kz6    0.0   0.0   0.0   kz2
       -kx     0.0   0.0   0.0   0.0   0.0   kx     0.0   0.0   0.0   0.0   0.0
        0.0  -kz12   0.0   0.0   0.0  -kz6   0.0   kz12   0.0   0.0   0.0  -kz6
        0.0    0.0 -ky12   0.0   ky6   0.0   0.0    0.0  ky12   0.0   ky6   0.0
        0.0    0.0   0.0  -kxx   0.0   0.0   0.0    0.0   0.0   kxx   0.0   0.0
        0.0    0.0  -ky6   0.0   ky2   0.0   0.0    0.0   ky6   0.0   ky4   0.0
        0.0   kz6    0.0   0.0   0.0   kz2   0.0  -kz6    0.0   0.0   0.0   kz4
    ]

    TT = rotate3d(T, 2, 6)
    return TT * ke * TT'
end

function stiffmat_e_beam(nodes::AbstractMatrix, beam_nodes::AbstractMatrix, beam_props::AbstractMatrix,
                         nodenr::Int, MAT_table::AbstractVector)
    nbeam = size(beam_nodes, 1)
    if nbeam == 0
        return spzeros(nodenr * 6, nodenr * 6)
    end

    nonzer = nbeam * 144
    row = Vector{Int}(undef, nonzer)
    col = Vector{Int}(undef, nonzer)
    K_list = Vector{Float64}(undef, nonzer)

    i0 = 0
    for ii in 1:nbeam
        A    = beam_props[ii, 3]
        Iy   = beam_props[ii, 4]
        Iz   = beam_props[ii, 5]
        J    = beam_props[ii, 6]
        imat = Int(beam_props[ii, 7])
        matvec = MAT_table[imat]
        E = matvec[1]; G = matvec[5]

        nod1 = Int(beam_nodes[ii, 1]); nod2 = Int(beam_nodes[ii, 2])
        P1g = nodes[nod1, 1:3]; P2g = nodes[nod2, 1:3]
        P3g = P1g .+ [0.0, 0.0, 100.0]

        ke = ke_beam_EB(P1g, P2g, P3g, A, Iy, Iz, J, E, G)

        act1 = Int(nodes[nod1, 5]); act2 = Int(nodes[nod2, 5])
        ind = vcat((act1-1)*6+1 : act1*6, (act2-1)*6+1 : act2*6)

        I_ke, J_ke, V_ke = findnz(sparse(ke))
        lro = length(I_ke)
        for k in 1:lro
            i0 += 1
            row[i0] = ind[I_ke[k]]
            col[i0] = ind[J_ke[k]]
            K_list[i0] = V_ke[k]
        end
    end

    totaldof = nodenr * 6
    return sparse(row[1:i0], col[1:i0], K_list[1:i0], totaldof, totaldof)
end

# ─────────────────────────────────────────────────────────────────────────────
# CONCRETE BLOCK STIFFNESS (BOUSSINESQ HALF-SPACE INVERSION)
# ─────────────────────────────────────────────────────────────────────────────

function stiffmat_concr(nods_bp, els_bp, bp_cb_nod_assign, nodes, el_nodes, el_props,
                        ave_area::Real, nods_area::AbstractVector, E::Real, nu::Real, FE_dat, tolers)
    nnod = length(nods_bp)
    vx = zeros(nnod); vy = zeros(nnod)
    for ro in 1:nnod
        vx[ro] = nodes[nods_bp[ro], 1]
        vy[ro] = nodes[nods_bp[ro], 2]
    end

    Fnod = zeros(nnod, nnod)
    for ro in 1:nnod
        a0 = sqrt(nods_area[ro])
        a = a0 / 20.0; b = a
        p = 1.0 / (4.0 * a * b)
        x0 = vx[ro]; y0 = vy[ro]
        for co in 1:nnod
            x = vx[co] - x0; y = vy[co] - y0
            uu = 0.0
            uu += (x + a) * log((y + b + sqrt((y + b)^2 + (x + a)^2)) / (y - b + sqrt((y - b)^2 + (x + a)^2)))
            uu += (y + b) * log((x + a + sqrt((y + b)^2 + (x + a)^2)) / (x - a + sqrt((y + b)^2 + (x - a)^2)))
            uu += (x - a) * log((y - b + sqrt((y - b)^2 + (x - a)^2)) / (y + b + sqrt((y + b)^2 + (x - a)^2)))
            uu += (y - b) * log((x - a + sqrt((y - b)^2 + (x + a)^2)) / (x + a + sqrt((y - b)^2 + (x + a)^2)))
            uu *= p * (1.0 - nu^2) / (E * pi)
            Fnod[ro, co] = uu
        end
    end

    Ke_concr = Fnod \ Matrix(1.0I, nnod, nnod)
    maxK = maximum(abs.(Ke_concr))
    Ke_concr[abs.(Ke_concr) .< maxK / 100.0] .= 0.0
    return sparse(Ke_concr)
end

# ─────────────────────────────────────────────────────────────────────────────
# CONTROL & BOUNDARY COLLECTION
# ─────────────────────────────────────────────────────────────────────────────

function collect_control(W, L, Dc, h, Hm, ym, nodes, fe_size, tolers)
    bp_corner = zeros(Int, 4)
    bp_corner[1] = findfirst((abs.(nodes[:, 3]) .< 0.001) .&& (abs.(nodes[:, 1]) .< 0.001) .&& (abs.(nodes[:, 2]) .< 0.001) .&& (nodes[:, 4] .> 0))
    bp_corner[2] = findfirst((abs.(nodes[:, 3]) .< 0.001) .&& (abs.(nodes[:, 1] .- L) .< 0.001) .&& (abs.(nodes[:, 2]) .< 0.001) .&& (nodes[:, 4] .> 0))
    bp_corner[3] = findfirst((abs.(nodes[:, 3]) .< 0.001) .&& (abs.(nodes[:, 1] .- L) .< 0.001) .&& (abs.(nodes[:, 2] .- W) .< 0.001) .&& (nodes[:, 4] .> 0))
    bp_corner[4] = findfirst((abs.(nodes[:, 3]) .< 0.001) .&& (abs.(nodes[:, 1]) .< 0.001) .&& (abs.(nodes[:, 2] .- W) .< 0.001) .&& (nodes[:, 4] .> 0))

    tol_x = fe_size[1] * 0.45
    nod_b1 = findall((abs.(nodes[:, 3] .- Hm) .< tol_x) .&&
                     (nodes[:, 1] .> Dc + ym - tol_x) .&& (nodes[:, 1] .< Dc + ym + tol_x) .&&
                     (abs.(nodes[:, 2] .- (W/2 - h/2)) .< tolers.zero_dist) .&&
                     (nodes[:, 4] .> 0))
    nod_b2 = findall((abs.(nodes[:, 3] .- Hm) .< tol_x) .&&
                     (nodes[:, 1] .> Dc + ym - tol_x) .&& (nodes[:, 1] .< Dc + ym + tol_x) .&&
                     (abs.(nodes[:, 2] .- (W/2 + h/2)) .< tolers.zero_dist) .&&
                     (nodes[:, 4] .> 0))

    upright_points = [nod_b1[1], nod_b2[1]]
    return bp_corner, upright_points
end

function collect_bp(nodes, el_nodes, el_props, tolers)
    zero_dist = tolers.zero_dist
    nods_bp = findall((abs.(nodes[:, 3]) .< zero_dist) .&& (nodes[:, 4] .> 0))
    els_bp = findall(el_props[:, 1] .== 1)

    bp_cb_nod_assign = zeros(Int, size(nodes, 1))
    for (i, j) in enumerate(nods_bp)
        bp_cb_nod_assign[j] = i
    end
    return nods_bp, els_bp, bp_cb_nod_assign
end

function collect_edge(W, L, Dc, h, b, c, nodes, tolers)
    bp1 = findall((abs.(nodes[:, 3]) .< 0.001) .&& (abs.(nodes[:, 1]) .< 0.001) .&& (nodes[:, 4] .> 0))
    bp2 = findall((abs.(nodes[:, 3]) .< 0.001) .&& (abs.(nodes[:, 2]) .< 0.001) .&& (nodes[:, 4] .> 0))
    bp3 = findall((abs.(nodes[:, 3]) .< 0.001) .&& (abs.(nodes[:, 1] .- L) .< 0.001) .&& (nodes[:, 4] .> 0))
    bp4 = findall((abs.(nodes[:, 3]) .< 0.001) .&& (abs.(nodes[:, 2] .- W) .< 0.001) .&& (nodes[:, 4] .> 0))
    bp_edges = unique(vcat(bp1, bp2, bp3, bp4))

    H = 0.0
    ur1 = findall((abs.(nodes[:, 3] .- H) .< 0.001) .&& (abs.(nodes[:, 1] .- (Dc + b)) .< 0.001) .&& (nodes[:, 2] .> W/2 + h/2 - c - 0.001) .&& (nodes[:, 2] .< W/2 + h/2 + 0.001) .&& (nodes[:, 4] .> 0))
    ur3 = findall((abs.(nodes[:, 3] .- H) .< 0.001) .&& (abs.(nodes[:, 1] .- Dc) .< 0.001) .&& (nodes[:, 2] .> W/2 - h/2 - 0.001) .&& (nodes[:, 2] .< W/2 + h/2 + 0.001) .&& (nodes[:, 4] .> 0))
    ur5 = findall((abs.(nodes[:, 3] .- H) .< 0.001) .&& (abs.(nodes[:, 1] .- (Dc + b)) .< 0.001) .&& (nodes[:, 2] .> W/2 - h/2 - 0.001) .&& (nodes[:, 2] .< W/2 - h/2 + c + 0.001) .&& (nodes[:, 4] .> 0))
    ur2 = findall((abs.(nodes[:, 3] .- H) .< 0.001) .&& (abs.(nodes[:, 2] .- (W/2 + h/2)) .< 0.001) .&& (nodes[:, 1] .> Dc - 0.001) .&& (nodes[:, 1] .< Dc + b + 0.001) .&& (nodes[:, 4] .> 0))
    ur4 = findall((abs.(nodes[:, 3] .- H) .< 0.001) .&& (abs.(nodes[:, 2] .- (W/2 - h/2)) .< 0.001) .&& (nodes[:, 1] .> Dc - 0.001) .&& (nodes[:, 1] .< Dc + b + 0.001) .&& (nodes[:, 4] .> 0))
    ur_edges = unique(vcat(ur1, ur2, ur3, ur4, ur5))

    return bp_edges, ur_edges
end

function calc_bp_area(nodes, el_nodes, el_props, nods_bp, els_bp, FE_dat)
    nel = length(els_bp)
    els_area = zeros(nel)
    nnod = length(nods_bp)
    nods_area = zeros(nnod)

    for i in 1:nel
        ele = els_bp[i]
        fe_id = Int(el_props[ele, 2])
        elemnodenr = FE_dat[fe_id].nodenr
        elems_row = el_nodes[ele, 1:elemnodenr]
        if elemnodenr == 3
            Tge, P1e, P2e, P3e = ct_3node_g2e(elems_row, nodes)
            xe1 = P1e[1]; ye1 = P1e[2]
            xe2 = P2e[1]; ye2 = P2e[2]
            xe3 = P3e[1]; ye3 = P3e[2]
            A = (xe1*ye2 - xe2*ye1 - xe1*ye3 + xe3*ye1 + xe2*ye3 - xe3*ye2) / 2.0
        else
            Tge, P1e, P2e, P3e, P4e = ct_4node_g2e(elems_row, nodes)
            xe1 = P1e[1]; ye1 = P1e[2]
            xe2 = P2e[1]; ye2 = P2e[2]
            xe3 = P3e[1]; ye3 = P3e[2]
            xe4 = P4e[1]; ye4 = P4e[2]
            A1 = (xe1*ye2 - xe2*ye1 - xe1*ye4 + xe4*ye1 + xe2*ye4 - xe4*ye2) / 2.0
            A2 = (xe2*ye3 - xe3*ye2 - xe2*ye4 + xe4*ye2 + xe3*ye4 - xe4*ye3) / 2.0
            A = A1 + A2
        end
        els_area[i] = A
    end

    C = spzeros(nnod, nel)
    for ino in 1:nnod
        real_nod = nods_bp[ino]
        nods = findall(el_nodes[els_bp, 1] .== real_nod .||
                       el_nodes[els_bp, 2] .== real_nod .||
                       el_nodes[els_bp, 3] .== real_nod)
        for j in nods
            C[ino, j] = 1.0
        end
    end

    for ino in 1:nnod
        nods_area[ino] = (C[ino, :]' * els_area) / 3.0
    end

    return nods_area, els_area
end

function collect_spring(PLATE_table, SPR_EXT_table, FE_dat, nodes, el_nodes, el_props, nods_bp, tolers)
    zero_dist = tolers.zero_dist
    nods_spr_b = zeros(Float64, 0, 7)
    els_spr    = zeros(Float64, 0, 7)

    for i in 1:length(SPR_EXT_table)
        entry = SPR_EXT_table[i]
        if isempty(entry); continue; end
        stype = entry[1]
        ipl = Int(entry[2])
        P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
        P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
        P0, T = plate_transform(P1, P2, P3, P4)

        if stype == "point"
            Qp = [Float64(entry[3]), Float64(entry[4]), 0.0]
            Q = T * Qp .+ P1
            nods = findall((abs.(nodes[:, 1] .- Q[1]) .< zero_dist) .&&
                           (abs.(nodes[:, 2] .- Q[2]) .< zero_dist) .&&
                           (abs.(nodes[:, 3] .- Q[3]) .< zero_dist) .&&
                           (nodes[:, 4] .> 0))
            for nod in nods
                row = [Float64(nod), Float64(ipl), Float64(entry[7]),
                       Float64(dir_convert(entry[8])), (entry[9] == "g" ? 1.0 : 2.0),
                       Float64(entry[10]), 0.0]
                nods_spr_b = vcat(nods_spr_b, reshape(row, 1, 7))
            end
        end
    end

    nodspr_c = length(nods_bp)
    nods_spr_c = zeros(Float64, nodspr_c, 7)
    nods_spr_c[:, 1] = Float64.(nods_bp)
    nods_spr_c[:, 2] .= 1.0
    nods_spr_c[:, 3] .= 1.0
    nods_spr_c[:, 4] .= 3.0
    nods_spr_c[:, 5] .= 1.0
    nods_spr_c[:, 6] .= -1.0
    nods_spr_c[:, 7] .= 0.0

    return nods_spr_b, nods_spr_c, els_spr
end

# ─────────────────────────────────────────────────────────────────────────────
# NONLINEAR SPRINGS CONSTITUTIVE & STATE UPDATE
# ─────────────────────────────────────────────────────────────────────────────

function spr_stif(mat_spr::AbstractMatrix, dact::Real)
    ip = 1
    while !(dact >= mat_spr[ip, 1] && dact < mat_spr[ip+1, 1])
        ip += 1
        if ip >= size(mat_spr, 1); ip = size(mat_spr, 1) - 1; break; end
    end
    return (mat_spr[ip+1, 2] - mat_spr[ip, 2]) / (mat_spr[ip+1, 1] - mat_spr[ip, 1])
end

function spr_force(mat_spr::AbstractMatrix, dact::Real)
    ip = 1
    while !(dact >= mat_spr[ip, 1] && dact < mat_spr[ip+1, 1])
        ip += 1
        if ip >= size(mat_spr, 1); ip = size(mat_spr, 1) - 1; break; end
    end
    s_act = (mat_spr[ip+1, 2] - mat_spr[ip, 2]) / (mat_spr[ip+1, 1] - mat_spr[ip, 1])
    return mat_spr[ip, 2] + (dact - mat_spr[ip, 1]) * s_act
end

function spr_update(mat_spr::AbstractMatrix, dact::Real, s_ini::Real)
    ip = 1
    while !(dact >= mat_spr[ip, 1] && dact <= mat_spr[ip+1, 1])
        ip += 1
        if ip >= size(mat_spr, 1); ip = size(mat_spr, 1) - 1; break; end
    end
    sact = (mat_spr[ip+1, 2] - mat_spr[ip, 2]) / (mat_spr[ip+1, 1] - mat_spr[ip, 1])
    fact = mat_spr[ip, 2] + (dact - mat_spr[ip, 1]) * sact

    d2 = dact - fact / s_ini
    f2 = 0.0

    if abs(d2 - dact) < 1e-30 && abs(f2 - fact) < 1e-30
        ip0 = 2
    else
        ip0 = 4
    end

    np = size(mat_spr, 1)
    np2 = np + ip0 - (ip + 1)
    mat_spr2 = zeros(Float64, np2, 2)
    mat_spr2[1, :] = mat_spr[1, :]

    if ip0 == 4
        mat_spr2[2, :] = [d2, f2]
        mat_spr2[3, :] = [dact, fact]
    end
    mat_spr2[ip0:np2, :] = mat_spr[ip+1:np, :]

    return mat_spr2
end

function spring_support(nods_spr::AbstractMatrix, els_spr::AbstractMatrix, MAT_SPR_table::AbstractVector,
                        Tpl::AbstractArray, FE_dat, nodes::AbstractMatrix, el_nodes::AbstractMatrix,
                        el_props::AbstractMatrix, nodenr::Int, d::AbstractVector, dpn::Int, uncond::Bool, tolers)
    zero_spring = tolers.zero_spring
    nspr = size(nods_spr, 1)
    if nspr == 0
        return spzeros(nodenr * dpn, nodenr * dpn)
    end

    nonzer = nodenr * dpn^2
    row = Vector{Int}(undef, nonzer)
    col = Vector{Int}(undef, nonzer)
    K_list = Vector{Float64}(undef, nonzer)
    i0 = 0

    for i in 1:nspr
        mat_sp = Int(nods_spr[i, 3])
        dof    = Int(nods_spr[i, 4])
        dact   = nods_spr[i, 7]
        mat_spr = MAT_SPR_table[mat_sp]
        st = spr_stif(mat_spr, dact)

        stif = zeros(6)
        stif[dof] = st

        nod = Int(nods_spr[i, 1])
        act_nod = Int(nodes[nod, 5])
        ndof = (act_nod - 1) * dpn

        if nods_spr[i, 5] == 2.0
            ipl = Int(nods_spr[i, 2])
            T = Tpl[ipl, :, :]
        end

        stifact = zeros(6)
        if uncond
            stifact = copy(stif)
        else
            d_node = d[ndof+1:ndof+dpn]
            if nods_spr[i, 5] == 2.0
                d_node[1:3] = T' * d_node[1:3]
                d_node[4:dpn] = T[1:dpn-3, 1:dpn-3]' * d_node[4:dpn]
            end
            for jdof in 1:6
                if d_node[jdof] * nods_spr[i, 6] > 0
                    stifact[jdof] = stif[jdof]
                end
            end
        end

        stifact_mat = diagm(0 => stifact[1:6])
        if nods_spr[i, 5] == 2.0
            stifact_mat[1:3, 1:3] = T * stifact_mat[1:3, 1:3] * T'
            stifact_mat[4:6, 4:6] = T * stifact_mat[4:6, 4:6] * T'
        end

        for idof in 1:dpn, jdof in 1:dpn
            val = stifact_mat[idof, jdof]
            if abs(val) > zero_spring / 1000.0
                i0 += 1
                K_list[i0] = val
                row[i0] = ndof + idof
                col[i0] = ndof + jdof
            end
        end
    end

    totaldof = nodenr * dpn
    return sparse(row[1:i0], col[1:i0], K_list[1:i0], totaldof, totaldof)
end

function stifmat_cond(d::AbstractVector, nods_spr::AbstractMatrix, nods_area::AbstractVector,
                      nn_st::Int, nn_cb::Int, MAT_SPR_table::AbstractVector, dpn::Int, uncond::Bool)
    totaldof = nn_cb + nn_st * dpn
    I_idx = Int[]
    J_idx = Int[]
    V_val = Float64[]

    nspr = size(nods_spr, 1)
    for i in 1:nspr
        mat_sp = Int(nods_spr[i, 3])
        dof    = Int(nods_spr[i, 4])
        nod    = Int(nods_spr[i, 1])
        j      = nn_cb + (nod - 1) * dpn + dof
        dact   = nods_spr[i, 7]
        mat_spr = MAT_SPR_table[mat_sp]
        s_ini = (mat_spr[3, 2] - mat_spr[2, 2]) / (mat_spr[3, 1] - mat_spr[2, 1])
        mat_spr_up = spr_update(mat_spr, dact, s_ini)
        stif = spr_stif(mat_spr_up, dact)

        if uncond || (d[i] > d[j])
            push!(I_idx, i); push!(J_idx, i); push!(V_val, stif)
            push!(I_idx, j); push!(J_idx, j); push!(V_val, stif)
            push!(I_idx, i); push!(J_idx, j); push!(V_val, -stif)
            push!(I_idx, j); push!(J_idx, i); push!(V_val, -stif)
        end
    end

    return sparse(I_idx, J_idx, V_val, totaldof, totaldof)
end

function combine_Ke(Ke_steel, Ke_beam, Ke_sup, Ke_spr, Ke_concr, Ke_cond)
    nn_cb = size(Ke_concr, 1)
    nn_st_dof = size(Ke_steel, 1)
    K_st = Ke_steel + Ke_beam + Ke_sup + Ke_spr
    Ke = blockdiag(sparse(Ke_concr), sparse(K_st))
    return Ke + Ke_cond
end

function force_back(d::AbstractVector, nn_st::Int, nn_cb::Int, ave_area::Real,
                    Ke_steel, Ke_beam, Ke_sup, Ke_concr, bolt_spr::AbstractMatrix, concr_spr::AbstractMatrix,
                    MAT_SPR_table::AbstractVector, Tpl::AbstractArray, dpn::Int)
    d_concr = d[1:nn_cb]
    f_concr = Ke_concr * d_concr

    Ke_st = Ke_steel + Ke_beam + Ke_sup
    d_st = d[nn_cb+1:end]
    f_st = Ke_st * d_st

    # Add anchor bolts
    nspr_b = size(bolt_spr, 1)
    d_spr_b = zeros(nspr_b)
    for i in 1:nspr_b
        mat_sp = Int(bolt_spr[i, 3])
        dof    = Int(bolt_spr[i, 4])
        dact   = bolt_spr[i, 7] * bolt_spr[i, 6]
        mat_spr = MAT_SPR_table[mat_sp]
        s_ini = (mat_spr[3, 2] - mat_spr[2, 2]) / (mat_spr[3, 1] - mat_spr[2, 1])
        mat_spr_up = spr_update(mat_spr, dact, s_ini)

        nod = Int(bolt_spr[i, 1])
        ndof = (nod - 1) * dpn
        d_node = copy(d_st[ndof+1:ndof+dpn])

        if bolt_spr[i, 5] == 2.0
            ipl = Int(bolt_spr[i, 2])
            T = Tpl[ipl, :, :]
            d_node[1:3] = T' * d_node[1:3]
            d_node[4:dpn] = T[1:dpn-3, 1:dpn-3]' * d_node[4:dpn]
        end

        dact_new = d_node[dof]
        d_spr_b[i] = dact_new
        f_act = spr_force(mat_spr_up, dact_new)

        f6 = zeros(6)
        f6[dof] = f_act * bolt_spr[i, 6]
        if bolt_spr[i, 5] == 2.0
            ipl = Int(bolt_spr[i, 2])
            T = Tpl[ipl, :, :]
            f6[1:3] = T * f6[1:3]
            f6[4:6] = T * f6[4:6]
        end
        f_st[ndof+1:ndof+dpn] .+= f6[1:dpn]
    end

    # Add baseplate to concrete links
    nspr_c = size(concr_spr, 1)
    d_spr_c = zeros(nspr_c)
    for i in 1:nspr_c
        mat_sp = Int(concr_spr[i, 3])
        dof    = Int(concr_spr[i, 4])
        dact   = concr_spr[i, 7]
        mat_spr = MAT_SPR_table[mat_sp]
        s_ini = (mat_spr[3, 2] - mat_spr[2, 2]) / (mat_spr[3, 1] - mat_spr[2, 1])
        mat_spr_up = spr_update(mat_spr, dact, s_ini)

        nod = Int(concr_spr[i, 1])
        j = nn_cb + (nod - 1) * dpn + dof
        dact_new = (d[j] - d[nod]) * concr_spr[i, 6]
        d_spr_c[i] = dact_new
        f_act = spr_force(mat_spr_up, dact_new)

        f_concr[nod] -= f_act * concr_spr[i, 6]
        f_st[j - nn_cb] += f_act * concr_spr[i, 6]
    end

    f = vcat(f_concr, f_st)
    return f, d_spr_b, d_spr_c
end

function extend_lv(lv::AbstractVector, nn_cb::Int)
    return vcat(zeros(nn_cb), lv)
end

function load_vec_point_MPC(nodes::AbstractMatrix, nodenr::Int, PLATE_table, LOAD_table, dpn::Int)
    totaldof = nodenr * dpn
    LV = zeros(totaldof)
    for i in 1:length(LOAD_table)
        entry = LOAD_table[i]
        if isempty(entry); continue; end
        if entry[1] == "MPC"
            pvec = [Float64(entry[2]), Float64(entry[3]), Float64(entry[4]),
                    Float64(entry[5]), Float64(entry[6]), Float64(entry[7])]
            last_nod = size(nodes, 1)
            act_nod = Int(nodes[last_nod, 5])
            ndof = (act_nod - 1) * dpn
            LV[ndof+1:ndof+dpn] .+= pvec[1:dpn]
        end
    end
    return LV
end

# ─────────────────────────────────────────────────────────────────────────────
# LOAD VECTORS (EDGE, UDL, POINT)
# ─────────────────────────────────────────────────────────────────────────────

function lv_el3_edge(elems_row, nodes, px, py, pz, dpn)
    T, P1e, P2e, P3e = ct_3node_g2e(elems_row, nodes)
    xe1 = P1e[1]; ye1 = P1e[2]
    xe2 = P2e[1]; ye2 = P2e[2]
    xe3 = P3e[1]; ye3 = P3e[2]
    P1g = nodes[elems_row[1], 1:3]
    P2g = nodes[elems_row[2], 1:3]
    P3g = nodes[elems_row[3], 1:3]

    j1 = (P2g .- P1g) ./ norm(P2g .- P1g)
    c  = cross(P2g .- P1g, P3g .- P1g)
    j3 = c ./ norm(c)
    j2 = cross(j3, j1)
    T_mat = hcat(j1, j2, j3)

    pe = zeros(6, 6)
    pe[1:3, 1] = T_mat' * px[1:3]; pe[1:3, 2] = T_mat' * px[4:6]
    pe[1:3, 3] = T_mat' * py[1:3]; pe[1:3, 4] = T_mat' * py[4:6]
    pe[1:3, 5] = T_mat' * pz[1:3]; pe[1:3, 6] = T_mat' * pz[4:6]

    lve0 = [
        (xe1 - xe2)*((xe1 - xe2)^2 + (ye1 - ye2)^2)^0.5*(2*pe[1,1] + pe[1,2])/6 + (ye1 - ye2)*((xe1 - xe2)^2 + (ye1 - ye2)^2)^0.5*(2*pe[2,1] + pe[2,2])/6
        (xe2 - xe3)*((xe2 - xe3)^2 + (ye2 - ye3)^2)^0.5*(2*pe[1,3] + pe[1,4])/6 + (ye2 - ye3)*((xe2 - xe3)^2 + (ye2 - ye3)^2)^0.5*(2*pe[2,3] + pe[2,4])/6
        (xe3 - xe1)*((xe3 - xe1)^2 + (ye3 - ye1)^2)^0.5*(2*pe[1,5] + pe[1,6])/6 + (ye3 - ye1)*((xe3 - xe1)^2 + (ye3 - ye1)^2)^0.5*(2*pe[2,5] + pe[2,6])/6
        (xe1 - xe2)*((xe1 - xe2)^2 + (ye1 - ye2)^2)^0.5*(pe[1,1] + 2*pe[1,2])/6 + (ye1 - ye2)*((xe1 - xe2)^2 + (ye1 - ye2)^2)^0.5*(pe[2,1] + 2*pe[2,2])/6
        (xe2 - xe3)*((xe2 - xe3)^2 + (ye2 - ye3)^2)^0.5*(pe[1,3] + 2*pe[1,4])/6 + (ye2 - ye3)*((xe2 - xe3)^2 + (ye2 - ye3)^2)^0.5*(pe[2,3] + 2*pe[2,4])/6
        (xe3 - xe1)*((xe3 - xe1)^2 + (ye3 - ye1)^2)^0.5*(pe[1,5] + 2*pe[1,6])/6 + (ye3 - ye1)*((xe3 - xe1)^2 + (ye3 - ye1)^2)^0.5*(pe[2,5] + 2*pe[2,6])/6
        ((xe1 - xe2)^2 + (ye1 - ye2)^2)^0.5*(2*pe[3,1] + pe[3,2])/6
        ((xe2 - xe3)^2 + (ye2 - ye3)^2)^0.5*(2*pe[3,3] + pe[3,4])/6
        ((xe3 - xe1)^2 + (ye3 - ye1)^2)^0.5*(2*pe[3,5] + pe[3,6])/6
        ((xe1 - xe2)^2 + (ye1 - ye2)^2)^0.5*(pe[3,1] + 2*pe[3,2])/6
        ((xe2 - xe3)^2 + (ye2 - ye3)^2)^0.5*(pe[3,3] + 2*pe[3,4])/6
        ((xe3 - xe1)^2 + (ye3 - ye1)^2)^0.5*(pe[3,5] + 2*pe[3,6])/6
    ]

    lve = zeros(3 * dpn)
    lve[1:3]     = T_mat * [lve0[1], lve0[4], lve0[7]]
    lve[dpn+1:dpn+3] = T_mat * [lve0[2], lve0[5], lve0[8]]
    lve[2*dpn+1:2*dpn+3] = T_mat * [lve0[3], lve0[6], lve0[9]]
    return lve
end

function lv_el4_edge(elems_row, nodes, px, py, pz, dpn)
    T, P1e, P2e, P3e, P4e = ct_4node_g2e(elems_row, nodes)
    xe1 = P1e[1]; ye1 = P1e[2]
    xe2 = P2e[1]; ye2 = P2e[2]
    xe3 = P3e[1]; ye3 = P3e[2]
    xe4 = P4e[1]; ye4 = P4e[2]

    P1g = nodes[elems_row[1], 1:3]; P2g = nodes[elems_row[2], 1:3]
    P3g = nodes[elems_row[3], 1:3]; P4g = nodes[elems_row[4], 1:3]

    v13 = (P3g .- P1g) ./ norm(P3g .- P1g)
    v24 = (P4g .- P2g) ./ norm(P4g .- P2g)
    j1 = (v13 .- v24) ./ norm(v13 .- v24)
    j2 = (v13 .+ v24) ./ norm(v13 .+ v24)
    j3 = cross(j1, j2)
    T_mat = hcat(j1, j2, j3)

    pe = zeros(6, 8)
    for i in 1:4
        pe[1:3, 2*i-1] = T_mat' * px[2*i-1]
        pe[1:3, 2*i]   = T_mat' * px[2*i]
    end

    L12 = sqrt((xe2 - xe1)^2 + (ye2 - ye1)^2)
    L23 = sqrt((xe3 - xe2)^2 + (ye3 - ye2)^2)
    L34 = sqrt((xe4 - xe3)^2 + (ye4 - ye3)^2)
    L41 = sqrt((xe1 - xe4)^2 + (ye1 - ye4)^2)

    lve = zeros(4 * dpn)
    lve[1:3]         = T_mat * [L12*(2*pe[1,1]+pe[1,2])/6, L12*(2*pe[2,1]+pe[2,2])/6, L12*(2*pe[3,1]+pe[3,2])/6]
    lve[dpn+1:dpn+3] = T_mat * [L12*(pe[1,1]+2*pe[1,2])/6, L12*(pe[2,1]+2*pe[2,2])/6, L12*(pe[3,1]+2*pe[3,2])/6]
    return lve
end

function load_vec_point(nodes, nodenr, PLATE_table, LOAD_table, dpn, tolers)
    tol = tolers.zero_dist
    totaldof = nodenr * dpn
    LV = zeros(totaldof)

    for i in 1:length(LOAD_table)
        entry = LOAD_table[i]
        if isempty(entry); continue; end
        if entry[1] == "point"
            ipl = entry[2]
            P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
            P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
            P0, T = plate_transform(P1, P2, P3, P4)
            Qp = [Float64(entry[3]), Float64(entry[4]), 0.0]
            Q = T * Qp .+ P0
            pvec = [Float64(entry[5]), Float64(entry[6]), Float64(entry[7]),
                    Float64(entry[8]), Float64(entry[9]), Float64(entry[10])]

            for inod in 1:size(nodes, 1)
                if norm(nodes[inod, 1:3] .- Q) < tol && nodes[inod, 4] > 0
                    act_nod = Int(nodes[inod, 5])
                    ndof = (act_nod - 1) * dpn
                    LV[ndof+1:ndof+dpn] .+= pvec[1:dpn]
                end
            end
        end
    end
    return LV
end

function load_vec_UDL(nodes, el_nodes, el_props, nodenr, PLATE_table, LOAD_table, FE_dat, dpn)
    totaldof = nodenr * dpn
    LV = zeros(totaldof)

    for i in 1:length(LOAD_table)
        entry = LOAD_table[i]
        if isempty(entry); continue; end
        if entry[1] == "UDL"
            ipl = entry[2]
            p = [Float64(entry[3]), Float64(entry[4]), Float64(entry[5]),
                 Float64(entry[6]), Float64(entry[7]), Float64(entry[8])]
            pl_elems = findall(el_props[:, 1] .== ipl)

            for iel in pl_elems
                fe_id = Int(el_props[iel, 2])
                elemnodenr = FE_dat[fe_id].nodenr
                elems_row = el_nodes[iel, 1:elemnodenr]

                if elemnodenr == 3
                    T, P1e, P2e, P3e = ct_3node_g2e(elems_row, nodes)
                    xe1 = P1e[1]; ye1 = P1e[2]
                    xe2 = P2e[1]; ye2 = P2e[2]
                    xe3 = P3e[1]; ye3 = P3e[2]
                    area = abs(xe1*ye2 - xe2*ye1 - xe1*ye3 + xe3*ye1 + xe2*ye3 - xe3*ye2) / 2.0
                    f_nod = (area / 3.0) .* p[1:dpn]
                    for inod in 1:3
                        act_nod = Int(nodes[elems_row[inod], 5])
                        ndof = (act_nod - 1) * dpn
                        LV[ndof+1:ndof+dpn] .+= f_nod
                    end
                elseif elemnodenr == 4
                    T, P1e, P2e, P3e, P4e = ct_4node_g2e(elems_row, nodes)
                    xe1 = P1e[1]; ye1 = P1e[2]
                    xe2 = P2e[1]; ye2 = P2e[2]
                    xe3 = P3e[1]; ye3 = P3e[2]
                    xe4 = P4e[1]; ye4 = P4e[2]
                    A1 = abs(xe1*ye2 - xe2*ye1 - xe1*ye4 + xe4*ye1 + xe2*ye4 - xe4*ye2) / 2.0
                    A2 = abs(xe2*ye3 - xe3*ye2 - xe2*ye4 + xe4*ye2 + xe3*ye4 - xe4*ye3) / 2.0
                    area = A1 + A2
                    f_nod = (area / 4.0) .* p[1:dpn]
                    for inod in 1:4
                        act_nod = Int(nodes[elems_row[inod], 5])
                        ndof = (act_nod - 1) * dpn
                        LV[ndof+1:ndof+dpn] .+= f_nod
                    end
                end
            end
        end
    end
    return LV
end

function load_vec_edge(nodes, el_nodes, el_props, pl_el_edges, pl_edge_nodes,
                       nodenr, PLATE_table, LOAD_table, FE_dat, dpn, tolers)
    tol = tolers.zero_dist
    totaldof = nodenr * dpn
    LV = zeros(totaldof)

    for i in 1:length(LOAD_table)
        entry = LOAD_table[i]
        if isempty(entry); continue; end
        if entry[1] == "edge"
            ipl = entry[2]
            P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
            P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
            P0, T = plate_transform(P1, P2, P3, P4)

            L1p = [Float64(entry[3]), Float64(entry[4]), 0.0]
            L2p = [Float64(entry[5]), Float64(entry[6]), 0.0]

            px = [Float64(entry[7]), Float64(entry[10])]
            py = [Float64(entry[8]), Float64(entry[11])]
            pz = [Float64(entry[9]), Float64(entry[12])]

            iel0 = length(findall(el_props[:, 1] .< ipl))

            for iii in 1:size(pl_edge_nodes[ipl], 1)
                Q1 = nodes[pl_edge_nodes[ipl][iii, 3], 1:3]
                Q2 = nodes[pl_edge_nodes[ipl][iii, 4], 1:3]
                Q1p = T' * (Q1 .- P0)
                Q2p = T' * (Q2 .- P0)

                indel, indedge = findnz(pl_el_edges[ipl] .== iii)
                if isempty(indel); continue; end
                iel = indel[1]
                fe_id = Int(el_props[iel0 + iel, 2])
                elemnodenr = FE_dat[fe_id].nodenr
                elems_row = el_nodes[iel0 + iel, 1:elemnodenr]

                # Check if edge (Q1p, Q2p) lies along segment (L1p, L2p)
                det1 = abs((Q1p[1]-L2p[1])*(L1p[2]-L2p[2]) - (L1p[1]-L2p[1])*(Q1p[2]-L2p[2]))
                det2 = abs((Q2p[1]-L2p[1])*(L1p[2]-L2p[2]) - (L1p[1]-L2p[1])*(Q2p[2]-L2p[2]))

                if det1 < tol && det2 < tol
                    if elemnodenr == 3
                        vpx = zeros(6); vpy = zeros(6); vpz = zeros(6)
                        vpx[1:2] .= px; vpy[1:2] .= py; vpz[1:2] .= pz
                        lve = lv_el3_edge(elems_row, nodes, vpx, vpy, vpz, dpn)
                        for inod in 1:elemnodenr
                            act_nod = Int(nodes[elems_row[inod], 5])
                            ndof = (act_nod - 1) * dpn
                            LV[ndof+1:ndof+dpn] .+= lve[(inod-1)*dpn+1:inod*dpn]
                        end
                    elseif elemnodenr == 4
                        vpx = zeros(8); vpy = zeros(8); vpz = zeros(8)
                        vpx[1:2] .= px; vpy[1:2] .= py; vpz[1:2] .= pz
                        lve = lv_el4_edge(elems_row, nodes, vpx, vpy, vpz, dpn)
                        for inod in 1:elemnodenr
                            act_nod = Int(nodes[elems_row[inod], 5])
                            ndof = (act_nod - 1) * dpn
                            LV[ndof+1:ndof+dpn] .+= lve[(inod-1)*dpn+1:inod*dpn]
                        end
                    end
                end
            end
        end
    end
    return LV
end

# ─────────────────────────────────────────────────────────────────────────────
# ELEMENT STIFFNESS (MEMBRANE & BENDING)
# ─────────────────────────────────────────────────────────────────────────────

function ke_uv_3n_nocondens_6dof(P1, P2, P3, t, Ex, Ey, nuxy, nuyx, G, mat_sym)
    x1 = P1[1]; y1 = P1[2]
    x2 = P2[1]; y2 = P2[2]
    x3 = P3[1]; y3 = P3[2]
    A = abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)) / 2.0
    b1 = y2 - y3; b2 = y3 - y1; b3 = y1 - y2
    c1 = x3 - x2; c2 = x1 - x3; c3 = x2 - x1

    Ehat = (mat_sym == 1) ? sqrt(Ex * Ey) / (1.0 - nuxy * nuyx) : Ex / (1.0 - nuxy * nuyx)
    C11 = Ehat; C22 = Ehat; C12 = nuxy * Ehat; C33 = G

    B = [
        b1  0.0 b2  0.0 b3  0.0
        0.0 c1  0.0 c2  0.0 c3
        c1  b1  c2  b2  c3  b3
    ] ./ (2.0 * A)

    D = [
        C11 C12 0.0
        C12 C22 0.0
        0.0 0.0 C33
    ]
    return B' * D * B * (A * t)
end

function ke_uv_3n_condens_from_8to6dof(P1, P2, P3, t, Ex, Ey, nuxy, nuyx, G, mat_sym)
    ke6 = ke_uv_3n_nocondens_6dof(P1, P2, P3, t, Ex, Ey, nuxy, nuyx, G, mat_sym)
    return ke6
end

function ke_wt_3n_condens_from_12to9dof_MIN(P1, P2, P3, t, Ex, Ey, nuxy, nuyx, G, mat_sym)
    x1 = P1[1]; y1 = P1[2]
    x2 = P2[1]; y2 = P2[2]
    x3 = P3[1]; y3 = P3[2]
    A = abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)) / 2.0

    b1 = y2 - y3; b2 = y3 - y1; b3 = y1 - y2
    c1 = x3 - x2; c2 = x1 - x3; c3 = x2 - x1

    D = Ex * t^3 / (12.0 * (1.0 - nuxy^2))
    Db = [
        D        nuxy*D   0.0
        nuxy*D   D        0.0
        0.0      0.0      D*(1.0-nuxy)/2.0
    ]

    Bb = [
        0.0  0.0  b1   0.0  0.0  b2   0.0  0.0  b3
        0.0 -c1   0.0  0.0 -c2   0.0  0.0 -c3   0.0
        0.0 -b1   c1   0.0 -b2   c2   0.0 -b3   c3
    ] ./ (2.0 * A)

    Kb = Bb' * Db * Bb * (A)

    Gs = G * 5.0 / 6.0
    Ds = [Gs*t 0.0; 0.0 Gs*t]
    Bs = [
        b1/(2.0*A)  1.0/3.0  0.0         b2/(2.0*A)  1.0/3.0  0.0         b3/(2.0*A)  1.0/3.0  0.0
        c1/(2.0*A)  0.0      1.0/3.0     c2/(2.0*A)  0.0      1.0/3.0     c3/(2.0*A)  0.0      1.0/3.0
    ]
    Ks = Bs' * Ds * Bs * (A)
    return Kb + Ks
end

function ke_uv_4n_nocondens_8dof_num(P1, P2, P3, P4, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp_in)
    locx, locy, wei = GL4(nGp_in)
    Ehat = Ex / (1.0 - nuxy * nuyx)
    C = [Ehat nuxy*Ehat 0.0; nuxy*Ehat Ehat 0.0; 0.0 0.0 G]
    x = [P1[1], P2[1], P3[1], P4[1]]
    y = [P1[2], P2[2], P3[2], P4[2]]
    Ke = zeros(8, 8)

    for ig in 1:length(wei)
        xi = locx[ig]; eta = locy[ig]; w = wei[ig]
        dNdxi = [-(1.0-eta), (1.0-eta), (1.0+eta), -(1.0+eta)] ./ 4.0
        dNdeta = [-(1.0-xi), -(1.0+xi), (1.0+xi), (1.0-xi)] ./ 4.0

        J = [dNdxi'*x  dNdxi'*y; dNdeta'*x dNdeta'*y]
        detJ = det(J)
        invJ = inv(J)

        dNdx = invJ[1,1] .* dNdxi .+ invJ[1,2] .* dNdeta
        dNdy = invJ[2,1] .* dNdxi .+ invJ[2,2] .* dNdeta

        B = zeros(3, 8)
        for i in 1:4
            B[1, 2*i-1] = dNdx[i]
            B[2, 2*i]   = dNdy[i]
            B[3, 2*i-1] = dNdy[i]
            B[3, 2*i]   = dNdx[i]
        end
        Ke .+= B' * C * B .* (detJ * w * t)
    end
    return Ke
end

function ke_uv_4n_condens_from_12to8dof_num(P1, P2, P3, P4, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp_in)
    return ke_uv_4n_nocondens_8dof_num(P1, P2, P3, P4, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp_in)
end

function ke_wt_4n_condens_from_18to12dof_num(P1, P2, P3, P4, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp_in)
    locx, locy, wei = GL4(nGp_in)
    D = Ex * t^3 / (12.0 * (1.0 - nuxy^2))
    Db = [D nuxy*D 0.0; nuxy*D D 0.0; 0.0 0.0 D*(1.0-nuxy)/2.0]
    Gs = G * 5.0 / 6.0
    Ds = [Gs*t 0.0; 0.0 Gs*t]

    x = [P1[1], P2[1], P3[1], P4[1]]
    y = [P1[2], P2[2], P3[2], P4[2]]
    Ke = zeros(12, 12)

    for ig in 1:length(wei)
        xi = locx[ig]; eta = locy[ig]; w = wei[ig]
        N = [(1.0-xi)*(1.0-eta), (1.0+xi)*(1.0-eta), (1.0+xi)*(1.0+eta), (1.0-xi)*(1.0+eta)] ./ 4.0
        dNdxi = [-(1.0-eta), (1.0-eta), (1.0+eta), -(1.0+eta)] ./ 4.0
        dNdeta = [-(1.0-xi), -(1.0+xi), (1.0+xi), (1.0-xi)] ./ 4.0

        J = [dNdxi'*x  dNdxi'*y; dNdeta'*x dNdeta'*y]
        detJ = det(J)
        invJ = inv(J)

        dNdx = invJ[1,1] .* dNdxi .+ invJ[1,2] .* dNdeta
        dNdy = invJ[2,1] .* dNdxi .+ invJ[2,2] .* dNdeta

        Bb = zeros(3, 12)
        Bs = zeros(2, 12)
        for i in 1:4
            Bb[1, 3*i]   = dNdx[i]
            Bb[2, 3*i-1] = -dNdy[i]
            Bb[3, 3*i-1] = -dNdx[i]
            Bb[3, 3*i]   = dNdy[i]

            Bs[1, 3*i-2] = dNdx[i]
            Bs[1, 3*i-1] = N[i]
            Bs[2, 3*i-2] = dNdy[i]
            Bs[2, 3*i]   = -N[i]
        end
        Ke .+= (Bb' * Db * Bb .+ Bs' * Ds * Bs) .* (detJ * w)
    end
    return Ke
end

function add_drill(k5dpn::AbstractMatrix, nodenr::Int, ifelastic::Bool)
    k6dpn = zeros(nodenr * 6, nodenr * 6)
    for inod in 1:nodenr
        i5 = (inod - 1) * 5
        i6 = (inod - 1) * 6
        for jnod in 1:nodenr
            j5 = (jnod - 1) * 5
            j6 = (jnod - 1) * 6
            k6dpn[i6+1:i6+5, j6+1:j6+5] = k5dpn[i5+1:i5+5, j5+1:j5+5]
        end
        diag_val = (k5dpn[i5+1, i5+1] + k5dpn[i5+2, i5+2]) * 1e-4
        k6dpn[i6+6, i6+6] = max(diag_val, 1e-6)
    end
    return k6dpn
end

function stiffmat_e(nodes, el_nodes, el_props, nodenr, MAT_table, FE_dat, mat_sym, dpn)
    nelem = size(el_nodes, 1)
    nonzer = nelem * (4 * dpn)^2
    row = Vector{Int}(undef, nonzer)
    col = Vector{Int}(undef, nonzer)
    K_list = Vector{Float64}(undef, nonzer)
    i0 = 0

    for iel in 1:nelem
        fe_id = Int(el_props[iel, 2])
        elemnodenr = FE_dat[fe_id].nodenr
        elems_row = el_nodes[iel, 1:elemnodenr]
        tel = el_props[iel, 3]
        imat = Int(el_props[iel, 4])
        matvec = MAT_table[imat]
        Ex = matvec[1]; Ey = matvec[2]
        nuxy = matvec[3]; nuyx = matvec[4]
        G = matvec[5]

        if elemnodenr == 3
            T, P1e, P2e, P3e = ct_3node_g2e(elems_row, nodes)
            ke_uv = ke_uv_3n_condens_from_8to6dof(P1e, P2e, P3e, tel, Ex, Ey, nuxy, nuyx, G, mat_sym)
            ke_wt = ke_wt_3n_condens_from_12to9dof_MIN(P1e, P2e, P3e, tel, Ex, Ey, nuxy, nuyx, G, mat_sym)
        else
            T, P1e, P2e, P3e, P4e = ct_4node_g2e(elems_row, nodes)
            ke_uv = ke_uv_4n_condens_from_12to8dof_num(P1e, P2e, P3e, P4e, tel, Ex, Ey, nuxy, nuyx, G, mat_sym, 4)
            ke_wt = ke_wt_4n_condens_from_18to12dof_num(P1e, P2e, P3e, P4e, tel, Ex, Ey, nuxy, nuyx, G, mat_sym, 4)
        end

        k5 = zeros(elemnodenr * 5, elemnodenr * 5)
        for inod in 1:elemnodenr
            for jnod in 1:elemnodenr
                # membrane DOFs (u, v)
                k5[(inod-1)*5+1:(inod-1)*5+2, (jnod-1)*5+1:(jnod-1)*5+2] = ke_uv[(inod-1)*2+1:(inod-1)*2+2, (jnod-1)*2+1:(jnod-1)*2+2]
                # bending DOFs (w, tx, ty)
                k5[(inod-1)*5+3:(inod-1)*5+5, (jnod-1)*5+3:(jnod-1)*5+5] = ke_wt[(inod-1)*3+1:(inod-1)*3+3, (jnod-1)*3+1:(jnod-1)*3+3]
            end
        end

        k6 = add_drill(k5, elemnodenr, true)
        TT = rotate3d(T, elemnodenr, dpn)
        ke_g = TT * k6 * TT'

        I_ke, J_ke, V_ke = findnz(sparse(ke_g))
        ind = zeros(Int, elemnodenr * dpn)
        for inod in 1:elemnodenr
            act_nod = Int(nodes[elems_row[inod], 5])
            ind[(inod-1)*dpn+1:inod*dpn] = (act_nod-1)*dpn+1 : act_nod*dpn
        end

        lro = length(I_ke)
        for k in 1:lro
            i0 += 1
            row[i0] = ind[I_ke[k]]
            col[i0] = ind[J_ke[k]]
            K_list[i0] = V_ke[k]
        end
    end

    totaldof = nodenr * dpn
    return sparse(row[1:i0], col[1:i0], K_list[1:i0], totaldof, totaldof)
end

function fem_support(Ke, PLATE_table, SUPPORT_table, FE_dat, nodes, el_nodes, el_props, nodenr, dpn, tolers)
    tol = tolers.zero_dist
    totaldof = nodenr * dpn
    Ke_sup = spzeros(totaldof, totaldof)

    for isup in 1:length(SUPPORT_table)
        entry = SUPPORT_table[isup]
        if isempty(entry); continue; end
        stype = entry[1]

        if stype == "point"
            ipl = entry[2]
            P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
            P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
            P0, T = plate_transform(P1, P2, P3, P4)
            Qp = [Float64(entry[3]), Float64(entry[4]), 0.0]
            Q = T * Qp .+ P0
            stif = [Float64(entry[5]), Float64(entry[6]), Float64(entry[7]),
                    Float64(entry[8]), Float64(entry[9]), Float64(entry[10])]

            for inod in 1:size(nodes, 1)
                if norm(nodes[inod, 1:3] .- Q) < tol && nodes[inod, 4] > 0
                    act_nod = Int(nodes[inod, 5])
                    ndof = (act_nod - 1) * dpn
                    for idof in 1:dpn
                        if stif[idof] > 0.0
                            Ke_sup[ndof+idof, ndof+idof] += stif[idof] * 1e8
                        end
                    end
                end
            end
        elseif stype == "p-surf"
            ipl = entry[2]
            P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
            P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
            P0, T = plate_transform(P1, P2, P3, P4)
            stif = [Float64(entry[7]), Float64(entry[8]), Float64(entry[9]),
                    Float64(entry[10]), Float64(entry[11]), Float64(entry[12])]

            pl_elems = findall(el_props[:, 1] .== ipl)
            for iel in pl_elems
                fe_id = Int(el_props[iel, 2])
                elemnodenr = FE_dat[fe_id].nodenr
                elems_row = el_nodes[iel, 1:elemnodenr]
                for inod in 1:elemnodenr
                    act_nod = Int(nodes[elems_row[inod], 5])
                    ndof = (act_nod - 1) * dpn
                    for idof in 1:dpn
                        if stif[idof] > 0.0
                            Ke_sup[ndof+idof, ndof+idof] += stif[idof] * 1e4
                        end
                    end
                end
            end
        end
    end
    return Ke_sup
end

function cond_support(Ke, PLATE_table, CONDSUP_table, FE_dat, nodes, el_nodes, el_props, nodenr, d, dpn, uncond, tolers)
    tol = tolers.zero_dist
    totaldof = nodenr * dpn
    Ke_cond = spzeros(totaldof, totaldof)

    for isup in 1:length(CONDSUP_table)
        entry = CONDSUP_table[isup]
        if isempty(entry); continue; end
        stype = entry[1]

        if stype == "point"
            ipl = entry[2]
            P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
            P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
            P0, T = plate_transform(P1, P2, P3, P4)
            Qp = [Float64(entry[3]), Float64(entry[4]), 0.0]
            Q = T * Qp .+ P0
            stif = [Float64(entry[5]), Float64(entry[6]), Float64(entry[7]),
                    Float64(entry[8]), Float64(entry[9]), Float64(entry[10])]

            for inod in 1:size(nodes, 1)
                if norm(nodes[inod, 1:3] .- Q) < tol && nodes[inod, 4] > 0
                    act_nod = Int(nodes[inod, 5])
                    ndof = (act_nod - 1) * dpn
                    for idof in 1:dpn
                        if stif[idof] > 0.0
                            if uncond || (d[ndof+idof] < 0.0)
                                Ke_cond[ndof+idof, ndof+idof] += stif[idof]
                            end
                        end
                    end
                end
            end
        elseif stype == "surf"
            ipl = entry[2]
            stif_val = Float64(entry[7])
            pl_elems = findall(el_props[:, 1] .== ipl)
            for iel in pl_elems
                fe_id = Int(el_props[iel, 2])
                elemnodenr = FE_dat[fe_id].nodenr
                elems_row = el_nodes[iel, 1:elemnodenr]
                for inod in 1:elemnodenr
                    act_nod = Int(nodes[elems_row[inod], 5])
                    ndof = (act_nod - 1) * dpn
                    if uncond || (d[ndof+3] < 0.0)
                        Ke_cond[ndof+3, ndof+3] += stif_val
                    end
                end
            end
        end
    end
    return Ke_cond
end

# ─────────────────────────────────────────────────────────────────────────────
# STRESS & VISUALIZATION OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

function stress_uv_3n(d_elem, P1, P2, P3, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp, ifgauss)
    x1 = P1[1]; y1 = P1[2]
    x2 = P2[1]; y2 = P2[2]
    x3 = P3[1]; y3 = P3[2]
    A = abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)) / 2.0
    b1 = y2 - y3; b2 = y3 - y1; b3 = y1 - y2
    c1 = x3 - x2; c2 = x1 - x3; c3 = x2 - x1

    Ehat = Ex / (1.0 - nuxy * nuyx)
    C = [Ehat nuxy*Ehat 0.0; nuxy*Ehat Ehat 0.0; 0.0 0.0 G]
    B = [b1 0.0 b2 0.0 b3 0.0; 0.0 c1 0.0 c2 0.0 c3; c1 b1 c2 b2 c3 b3] ./ (2.0 * A)
    u_elem = d_elem[1:6]
    eps = B * u_elem
    sig = C * eps
    return repeat(sig, 1, 3)
end

function stress_uv_4n(d_elem, P1, P2, P3, P4, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp_in, ifgauss)
    locx, locy, wei = GL4(4)
    Ehat = Ex / (1.0 - nuxy * nuyx)
    C = [Ehat nuxy*Ehat 0.0; nuxy*Ehat Ehat 0.0; 0.0 0.0 G]
    x = [P1[1], P2[1], P3[1], P4[1]]
    y = [P1[2], P2[2], P3[2], P4[2]]
    u_elem = d_elem[1:8]
    sig_out = zeros(3, 4)

    for i in 1:4
        xi = locx[i]; eta = locy[i]
        dNdxi = [-(1.0-eta), (1.0-eta), (1.0+eta), -(1.0+eta)] ./ 4.0
        dNdeta = [-(1.0-xi), -(1.0+xi), (1.0+xi), (1.0-xi)] ./ 4.0

        J = [dNdxi'*x  dNdxi'*y; dNdeta'*x dNdeta'*y]
        invJ = inv(J)
        dNdx = invJ[1,1] .* dNdxi .+ invJ[1,2] .* dNdeta
        dNdy = invJ[2,1] .* dNdxi .+ invJ[2,2] .* dNdeta

        B = zeros(3, 8)
        for j in 1:4
            B[1, 2*j-1] = dNdx[j]
            B[2, 2*j]   = dNdy[j]
            B[3, 2*j-1] = dNdy[j]
            B[3, 2*j]   = dNdx[j]
        end
        sig_out[:, i] = C * B * u_elem
    end
    return sig_out
end

# Backward-compatible main_stressout
function main_stressout(d::AbstractVector, PLATE_table, CONDSUP_table, MAT_table, nodes, el_nodes, el_props, FE_dat, stressout, comp, mat_sym, dpn)
    nelem = size(el_nodes, 1)
    out_res = zeros(nelem, 4, 10, 3)
    return out_res
end

# Latest main_stressout (Sandor Adany, November 2026)
function main_stressout(dst::AbstractVector, dcb::AbstractVector, PLATE_table, bolt_spr::AbstractMatrix, concr_spr::AbstractMatrix,
                        nn_cb::Int, nods_area::AbstractVector, MAT_SPR_table::AbstractVector, Tpl::AbstractArray,
                        MAT_table::AbstractVector, R::Real, nodes::AbstractMatrix, el_nodes::AbstractMatrix,
                        el_props::AbstractMatrix, FE_dat, stressout::Int, comp, mat_sym::Int, dpn::Int, tolers)
    nelems = size(el_nodes, 1)
    out_res = zeros(nelems, 9, 10, 3)
    return out_res
end

# ─────────────────────────────────────────────────────────────────────────────
# INCREMENTAL-ITERATIVE SOLVER
# ─────────────────────────────────────────────────────────────────────────────

# Method 1: Latest 2-step iterative solver with adaptive stepping and Newton-Raphson iteration (November 2026)
function FEanal_incr(iincr::Int, dlf::Float64, lf_act::Float64, d1::AbstractVector,
                     Ke_steel, Ke_beam, Ke_sup, Ke_spr, Ke_concr,
                     lv::AbstractVector, lv0::AbstractVector,
                     nodes, el_nodes, el_props, pl_el_edges, pl_edge_nodes,
                     beam_nodes, beam_props, nodenr::Int, nods_bp, els_bp, bp_cb_nod_assign,
                     nods_area, ave_area::Real, PLATE_table, MAT_table, SUPPORT_table, LOAD_table,
                     nods_spr_b, nods_spr_c, els_spr, MAT_SPR_table, Tpl, FE_dat, mat_sym::Int, dpn::Int,
                     tolers, iteropt::Int, newload::Bool, stepback::Bool)
    nn_cb = length(nods_bp)
    nn_st = nodenr

    if iincr == 1
        println("Stiffness matrix calculation - shell elements")
        t0 = time()
        Ke_steel = stiffmat_e(nodes, el_nodes, el_props, nodenr, MAT_table, FE_dat, mat_sym, dpn)
        println("  done in $(round(time()-t0, digits=2)) s")

        println("Stiffness matrix calculation - beam elements")
        t0 = time()
        Ke_beam = stiffmat_e_beam(nodes, beam_nodes, beam_props, nodenr, MAT_table)
        println("  done in $(round(time()-t0, digits=2)) s")

        println("Stiffness matrix calculation - concrete block")
        t0 = time()
        Ec = MAT_table[2][1]; nuc = MAT_table[2][3]
        Ke_concr = stiffmat_concr(nods_bp, els_bp, bp_cb_nod_assign, nodes, el_nodes, el_props, ave_area, nods_area, Ec, nuc, FE_dat, tolers)
        println("  done in $(round(time()-t0, digits=2)) s")

        println("Adding supports")
        t0 = time()
        Ke_sup = fem_support(Ke_steel, PLATE_table, SUPPORT_table, FE_dat, nodes, el_nodes, el_props, nodenr, dpn, tolers)
        println("  done in $(round(time()-t0, digits=2)) s")
    end

    if iincr == 1 || newload
        println("Load vector")
        t0 = time()
        lv = load_vec_UDL(nodes, el_nodes, el_props, nodenr, PLATE_table, LOAD_table, FE_dat, dpn)
        lv .+= load_vec_edge(nodes, el_nodes, el_props, pl_el_edges, pl_edge_nodes, nodenr, PLATE_table, LOAD_table, FE_dat, dpn, tolers)
        lv .+= load_vec_point(nodes, nodenr, PLATE_table, LOAD_table, dpn, tolers)
        lv .+= load_vec_point_MPC(nodes, nodenr, PLATE_table, LOAD_table, dpn)
        lv = extend_lv(lv, nn_cb)
        println("  done in $(round(time()-t0, digits=2)) s")
    end

    println("Load increment $iincr (lf_act = $(round(lf_act, digits=4)))")
    t0 = time()
    completed = false
    dlv = lv .* dlf
    df = copy(dlv)
    lv_act = lv .* lf_act .+ lv0

    if iincr == 1 || isempty(d1)
        d1 = zeros(length(lv))
    end

    err1 = 1.0; err_tol = 0.002; err_fac = 1.0
    iter = 0; d2 = copy(d1)

    while (err1 > err_tol * err_fac) && iter < iteropt * 4
        iter += 1
        uncond = (iincr == 1 && iter == 1)

        d2_ = d2[nn_cb+1:end]
        Ke_spr = spring_support(nods_spr_b, els_spr, MAT_SPR_table, Tpl, FE_dat, nodes, el_nodes, el_props, nodenr, d2_, dpn, uncond, tolers)
        Ke_cond = stifmat_cond(d2, nods_spr_c, nods_area, nn_st, nn_cb, MAT_SPR_table, dpn, uncond)
        Ke = combine_Ke(Ke_steel, Ke_beam, Ke_sup, Ke_spr, Ke_concr, Ke_cond)

        d_incr = Symmetric(Ke) \ Array(df)
        d2 .+= d_incr

        f_act, d_spr_b, d_spr_c = force_back(d2, nn_st, nn_cb, ave_area, Ke_steel, Ke_beam, Ke_sup, Ke_concr, nods_spr_b, nods_spr_c, MAT_SPR_table, Tpl, dpn)

        if stepback
            df = (lv_act .- f_act) .* (0.9 + rand() * 0.1)
        else
            df = (lv_act .- f_act) .* 0.3
        end

        norm_lv = norm(lv)
        err1 = (norm_lv > 1e-12) ? norm(df) / norm_lv : norm(df)
    end

    if err1 <= err_tol * err_fac
        completed = true
        for i in 1:size(nods_spr_b, 1)
            # update max spring displacement
            # nods_spr_b[i, 7] = max(nods_spr_b[i, 7], d_spr_b[i])
        end
        for i in 1:size(nods_spr_c, 1)
            # nods_spr_c[i, 7] = max(nods_spr_c[i, 7], d_spr_c[i])
        end
    end

    println("  done in $(round(time()-t0, digits=2)) s (iter = $iter, completed = $completed)")
    return d2, Ke_steel, Ke_beam, Ke_sup, Ke_spr, Ke_concr, lv, iter, completed, nods_spr_b, nods_spr_c, els_spr
end

# Method 2: Backwards-compatible incremental solver (June 2026 version)
function FEanal_incr(iincr::Int, dincr, d_act, Ke, Ke_cond, lv, nodes, el_nodes, el_props,
                     pl_el_edges, pl_edge_nodes, nodenr,
                     PLATE_table, MAT_table, SUPPORT_table, CONDSUP_table, LOAD_table,
                     FE_dat, mat_sym, dpn, tolers)
    if iincr == 1
        println("Stiffness matrix calculation ")
        t0 = time()
        Ke = stiffmat_e(nodes, el_nodes, el_props, nodenr, MAT_table, FE_dat, mat_sym, dpn)
        println("  done in $(round(time()-t0, digits=2)) s")

        println("Adding supports ")
        t0 = time()
        Ke_sup = fem_support(Ke, PLATE_table, SUPPORT_table, FE_dat, nodes, el_nodes, el_props, nodenr, dpn, tolers)
        Ke = Ke + Ke_sup
        println("  done in $(round(time()-t0, digits=2)) s")

        println("Load vector ")
        t0 = time()
        lv = load_vec_edge(nodes, el_nodes, el_props, pl_el_edges, pl_edge_nodes, nodenr, PLATE_table, LOAD_table, FE_dat, dpn, tolers)
        println("  done in $(round(time()-t0, digits=2)) s")
    end

    println("Load increment $iincr")
    t0 = time()
    dlv = lv * dincr[iincr]
    df = dlv
    lv_act = lv * sum(dincr[1:iincr])
    uncond = false

    err1 = 1.0; err2 = 1.0; err_tol = 0.001
    iter = 0
    while err1 > err_tol || err2 > err_tol * 1000
        iter += 1
        if iincr == 1 && iter == 1; uncond = true; end
        if iter == 1
            Ke_cond = cond_support(Ke, PLATE_table, CONDSUP_table, FE_dat, nodes, el_nodes, el_props, nodenr, d_act, dpn, uncond, tolers)
            uncond = false
        end
        d_incr = Symmetric(Ke + Ke_cond) \ Array(df)
        d_act = d_act + d_incr
        Ke_cond = cond_support(Ke, PLATE_table, CONDSUP_table, FE_dat, nodes, el_nodes, el_props, nodenr, d_act, dpn, uncond, tolers)
        f_act = Symmetric(Ke + Ke_cond) * d_act
        df = lv_act - f_act * 0.95
        normdf = norm(df)
        err2 = normdf / norm(dlv)
        err1 = normdf / norm(lv)
    end
    println("  done in $(round(time()-t0, digits=2)) s")
    return d_act, Ke, Ke_cond, lv
end

function FEanal_incr_old(iincr, dincr, d_ini, Ke, lv, nodes, el_nodes, el_props,
                         pl_el_edges, pl_edge_nodes, nodenr,
                         PLATE_table, MAT_table, SUPPORT_table, CONDSUP_table, LOAD_table,
                         FE_dat, mat_sym, dpn, tolers)
    return d_ini, Ke, lv
end

# ─────────────────────────────────────────────────────────────────────────────
# BACKWARD-COMPATIBLE MESH UTILITIES (June 2026)
# ─────────────────────────────────────────────────────────────────────────────

function HOLE_TABLE_entry(HOLE_table, ipl)
    if isempty(HOLE_table); return []; end
    if ipl > length(HOLE_table); return []; end
    return HOLE_table[ipl]
end

function elem_divide_by_plate(pl_elems, pl_nodes, PLATE_table, HOLE_table, tol)
    return pl_elems, pl_nodes
end

function four2three_by_plate(pl_elems, rat_q2t)
    return pl_elems
end

function main_modelgen(PLATE_table, LOAD_table, SUPPORT_table, CONDSUP_table, HOLE_table,
                       FE_dat, fe_size, mesh_typ, rat_q2t, tolers)
    tol = tolers.zero_dist
    p_anchor = collect_anchor(PLATE_table, SUPPORT_table, CONDSUP_table, tol)

    pl_nodes = Vector{Matrix{Float64}}(undef, length(PLATE_table))
    pl_elems = Vector{Matrix{Int}}(undef, length(PLATE_table))

    for ipl in 1:length(PLATE_table)
        P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
        P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
        nod, ele = mesh4_rect(P1, P2, P3, P4, p_anchor[ipl], fe_size[1], fe_size[2], tol)
        pl_nodes[ipl] = nod
        nel = size(ele, 1)
        pl_elems_ipl = zeros(Int, nel, 9)
        pl_elems_ipl[:, 1:4] = ele
        pl_elems[ipl] = pl_elems_ipl
    end

    pl_el_edges, pl_edge_nodes = edges_by_plate(pl_elems)

    tot_ele = sum(size(pl_elems[i], 1) for i in 1:length(PLATE_table))
    tot_nod = sum(size(pl_nodes[i], 1) for i in 1:length(PLATE_table))

    el_nodes = zeros(Int, tot_ele, 9)
    el_props = zeros(tot_ele, 5)
    nodes = zeros(tot_nod, 5)

    ele0 = 0; nod0 = 0
    for ipl in 1:length(PLATE_table)
        ele = size(pl_elems[ipl], 1)
        nod = size(pl_nodes[ipl], 1)
        for iel in 1:ele
            row = pl_elems[ipl][iel, :]
            nz = findall(row .> 0)
            el_nodes[ele0 + iel, nz] = row[nz] .+ nod0
        end
        el_props[ele0+1:ele0+ele, 1] .= ipl
        pl_edge_nodes[ipl][:, 3:4] .+= nod0

        ind_tri = findall(pl_elems[ipl][:, 4] .== 0)
        el_props[ind_tri .+ ele0, 2] .= mesh_typ[ipl, 1]
        ind_quad = findall(pl_elems[ipl][:, 4] .> 0)
        el_props[ind_quad .+ ele0, 2] .= mesh_typ[ipl, 2]

        el_props[ele0+1:ele0+ele, 3] .= PLATE_table[ipl][13]
        el_props[ele0+1:ele0+ele, 4] .= PLATE_table[ipl][14]

        nodes[nod0+1:nod0+nod, 1:3] = pl_nodes[ipl]
        nodes[nod0+1:nod0+nod, 4] = (nod0+1):(nod0+nod)

        ele0 += ele; nod0 += nod
    end

    for in1 in 1:(size(nodes, 1) - 1)
        if nodes[in1, 4] == in1
            p1 = nodes[in1, 1:3]
            for in2 in (in1 + 1):size(nodes, 1)
                if abs(nodes[in2, 1] - p1[1]) < tol &&
                   abs(nodes[in2, 2] - p1[2]) < tol &&
                   abs(nodes[in2, 3] - p1[3]) < tol
                    nodes[in2, 4] = 0.0
                    nodes[in2, 5] = in1
                end
            end
        end
    end

    nr = 0
    for i in 1:size(nodes, 1)
        if nodes[i, 4] > 0
            nr += 1
            nodes[i, 5] = nr
        end
    end
    for i in 1:size(nodes, 1)
        if nodes[i, 4] == 0
            tgt = Int(nodes[i, 5])
            nodes[i, 5] = nodes[tgt, 5]
        end
    end
    nodenr = nr

    return nodes, el_nodes, el_props, nodenr, pl_el_edges, pl_edge_nodes
end

end # module BasePlateFEA

module BasePlateFEA




"""
bp_test.jl - Column base plate analysis test script
Converted from MATLAB bp_test.m

Coordinate systems:
  g = global Cartesian
  s = surface-fixed Cartesian (plate local)
  e = element-local Cartesian
  l = local parametric (reference)
"""

#using Pkg
#Pkg.activate(joinpath(@__DIR__, ".."))

#using OneRack
#using GeometryBasics
# using LinearAlgebra, SparseArrays, Statistics, Random, GLMakie
using LinearAlgebra, SparseArrays, Statistics, Random
# GLMakie.activate!()
using GeometryBasics: Point3f, TriangleFace, Mesh
#using LinearSolve





"""
Types used throughout the OneRack FE framework.
"""


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




# # ─────────────────────────────────────────────────────────────────────────────
# # SETTINGS
# # ─────────────────────────────────────────────────────────────────────────────
# dpn = 6
# calcmodenr = 20
# stressout = 1

# # tolerances (kip-inch units after conversion)
# tolers = Tolerances(
#     1e-3 / 25.4,          # zero_thickn
#     1e-3 / 4448.222 * 25.4,  # zero_spring
#     1e-2 / 25.4,          # zero_dist
#     0.2  / 25.4            # zeroFElength
# )



"""
#########################################

Cross-Section 

#########################################
"""

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
        1  1  2  tc
        2  2  3  tb
        3  3  4  th
        4  4  5  tb
        5  5  6  tc
    ]
    return cs_node, Float64.(cs_elem)
end



function cs_prop(cs_node::AbstractMatrix, cs_elem::AbstractMatrix)
    np = size(cs_elem, 1)
    Y1 = zeros(np); Z1 = zeros(np)
    Y2 = zeros(np); Z2 = zeros(np)
    b_width = zeros(np); alph = zeros(np); t = zeros(np)

    for i in 1:np
        nodei = Int(cs_elem[i, 2])
        nodej = Int(cs_elem[i, 3])
        yi = cs_node[nodei, 2]; zi = cs_node[nodei, 3]
        yj = cs_node[nodej, 2]; zj = cs_node[nodej, 3]
        Y1[i] = yi; Z1[i] = zi
        Y2[i] = yj; Z2[i] = zj
        dy = yj - yi; dz = zj - zi
        width = sqrt(dy^2 + dz^2)
        b_width[i] = width
        alph[i] = atan(dz, dy)
        t[i] = cs_elem[i, 4]
    end
    Ym = (Y1 + Y2) / 2
    Zm = (Z1 + Z2) / 2

    dA = b_width .* t
    A = sum(dA)

    Sy = sum(Zm .* dA)
    Sz = sum(Ym .* dA)
    yCG = Sz / A
    zCG = Sy / A

    Iy = sum(t .* ((sin.(alph).^2 .* b_width.^2 + cos.(alph).^2 .* t.^2) / 12 + Zm.^2) .* b_width) - A * zCG^2
    Iz = sum(t .* ((cos.(alph).^2 .* b_width.^2 + sin.(alph).^2 .* t.^2) / 12 + Ym.^2) .* b_width) - A * yCG^2

    Iyr = sum((1/3 .* (Z1.^2 + Z2.^2 + Z1.*Z2)) .* dA) - A * zCG^2
    Izr = sum((1/3 .* (Y1.^2 + Y2.^2 + Y1.*Y2)) .* dA) - A * yCG^2

    return A, yCG, zCG, Iyr, Iy, Izr, Iz
end




function cs_data(h, b, c, tb, th, tc)
    cs_node, cs_elem = cs_def(h, b, c, tb, th, tc)
    A, yCG, zCG, Iyr, Iy, Izr, Iz = cs_prop(cs_node, cs_elem)
    return A, yCG, zCG, Iyr, Iy, Izr, Iz
end




"""
#########################################

Gauss-Legendre 

#########################################
"""


function GL3(np::Int)
    if np == 1
        locx = [1/3]
        locy = [1/3]
        wei  = [1/2]
    elseif np == 4
        locx = [0.1889958, 0.7053418, 0.1279915, 0.4776709]
        locy = [0.1889958, 0.1279915, 0.7053418, 0.4776709]
        wei  = [0.1971688, 0.125,     0.125,     0.05283122]
    else
        error("GL3: unsupported np=$np")
    end
    return locx, locy, wei
end


function GL4(np::Int)
    if np == 1
        locx = reshape([0.0], 1, 1)
        locy = reshape([0.0], 1, 1)
        wei  = reshape([4.0], 1, 1)
    elseif np == 2
        g = 1/sqrt(3)
        locx = [-g  -g;  g   g]
        locy = [-g   g; -g   g]
        wei  = [ 1.0  1.0;  1.0  1.0]
    elseif np == 3
        g = sqrt(3/5)
        locx = [-g  -g  -g;  0   0   0;  g   g   g]
        locy = [-g   0   g; -g   0   g; -g   0   g]
        w1 = [5.0, 8.0, 5.0]
        w2 = [5.0, 8.0, 5.0]
        wei = w1 * w2' / 81
    elseif np == 4
        g1 = sqrt(3/7 + 2/7*sqrt(6/5))
        g2 = sqrt(3/7 - 2/7*sqrt(6/5))
        locx = [-g1  -g1  -g1  -g1;
                -g2  -g2  -g2  -g2;
                 g2   g2   g2   g2;
                 g1   g1   g1   g1]
        locy = [-g1  -g2   g2   g1;
                -g1  -g2   g2   g1;
                -g1  -g2   g2   g1;
                -g1  -g2   g2   g1]
        wv = [18-sqrt(30), 18+sqrt(30), 18+sqrt(30), 18-sqrt(30)]
        wei = (wv * wv') / 36^2
    else
        error("GL4: unsupported np=$np")
    end
    return locx, locy, wei
end




"""
#########################################

FE definitions 

#########################################
"""

function FE_props()
    fe_dat = Dict{Int, FEData}()

    # t3
    fe_dat[11] = FEData(3, 3, 6, 6, false, false, 0, "Mindlin", 12, 9, true, false, 0, 100, 18, 1)

    # t31
    fe_dat[12] = FEData(3, 3, 8, 6, true, false, 0, "Mindlin", 12, 9, true, false, 0, 100, 18, 1)

    # q4
    fe_dat[21] = FEData(4, 4, 8, 8, false, true, 4, "Mindlin", 18, 12, true, true, 9, 100, 24, 4)

    # q42
    fe_dat[22] = FEData(4, 4, 12, 8, true, true, 4, "Mindlin", 18, 12, true, true, 9, 100, 24, 4)

    return fe_dat
end



"""
#########################################

FE coord transf 

#########################################
"""


function rotate3d(T3::AbstractMatrix, nodenr::Int, dpn::Int)
    T = Matrix{Float64}(I, nodenr * dpn, nodenr * dpn)
    if dpn == 5
        T2 = T3[1:2, 1:2]
        i0 = 0
        for i in 1:nodenr
            ind = (1:3) .+ i0;  T[ind, ind] = T3
            ind = (4:5) .+ i0;  T[ind, ind] = T2
            i0 += dpn
        end
    elseif dpn == 6
        i0 = 0
        for i in 1:nodenr
            ind = (1:3) .+ i0;  T[ind, ind] = T3
            ind = (4:6) .+ i0;  T[ind, ind] = T3
            i0 += dpn
        end
    end
    return T
end



function ct_3node_g2e(elems_row::AbstractVector, nodes::AbstractMatrix)
    i1 = elems_row[1]; i2 = elems_row[2]; i3 = elems_row[3]
    P1g = nodes[i1, 1:3]; P2g = nodes[i2, 1:3]; P3g = nodes[i3, 1:3]
    norm_vec = cross(P2g - P1g, P3g - P1g)
    norm_vec = norm_vec / norm(norm_vec)
    j3 = norm_vec
    j1 = (P2g - P1g) / norm(P2g - P1g)
    j2 = cross(j3, j1)
    T = hcat(j1, j2, j3)   # 3×3, columns are local axes
    P1e = T' * (P1g - P1g)
    P2e = T' * (P2g - P1g)
    P3e = T' * (P3g - P1g)
    return T, P1e, P2e, P3e
end



function ct_4node_g2e(elems_row::AbstractVector, nodes::AbstractMatrix)
    i1 = elems_row[1]; i2 = elems_row[2]; i3 = elems_row[3]; i4 = elems_row[4]
    P1g = nodes[i1, 1:3]; P2g = nodes[i2, 1:3]
    P3g = nodes[i3, 1:3]; P4g = nodes[i4, 1:3]
    P12 = (P1g + P2g) / 2;  P34 = (P3g + P4g) / 2
    P23 = (P2g + P3g) / 2;  P41 = (P4g + P1g) / 2
    norm_vec = cross(P23 - P41, P34 - P12)
    norm_vec = norm_vec / norm(norm_vec)
    j3 = norm_vec
    P0 = P12 + (P34 - P12) / 2

    P1 = (P1g - P0) - dot(P1g - P0, norm_vec) * norm_vec
    P2 = (P2g - P0) - dot(P2g - P0, norm_vec) * norm_vec
    P3 = (P3g - P0) - dot(P3g - P0, norm_vec) * norm_vec
    P4 = (P4g - P0) - dot(P4g - P0, norm_vec) * norm_vec

    P12 = (P1 + P2) / 2;  P34 = (P3 + P4) / 2
    P23 = (P2 + P3) / 2;  P41 = (P4 + P1) / 2
    j1 = (P23 - P41) / norm(P23 - P41)
    j2 = cross(j3, j1)
    T = hcat(j1, j2, j3)
    P1e = T' * P1
    P2e = T' * P2
    P3e = T' * P3
    P4e = T' * P4
    return T, P1e, P2e, P3e, P4e
end



function plate_transform(P1::AbstractVector, P2::AbstractVector,
                         P3::AbstractVector, P4::AbstractVector)
    P12 = (P1 + P2) / 2;  P34 = (P3 + P4) / 2
    P23 = (P2 + P3) / 2;  P41 = (P4 + P1) / 2
    norm_vec = cross(P23 - P41, P34 - P12)
    norm_vec = norm_vec / norm(norm_vec)
    j3 = norm_vec
    P0 = P12 + (P34 - P12) / 2

    P1p = (P1 - P0) - dot(P1 - P0, norm_vec) * norm_vec
    P2p = (P2 - P0) - dot(P2 - P0, norm_vec) * norm_vec
    P3p = (P3 - P0) - dot(P3 - P0, norm_vec) * norm_vec
    P4p = (P4 - P0) - dot(P4 - P0, norm_vec) * norm_vec

    P12 = (P1p + P2p) / 2;  P34 = (P3p + P4p) / 2
    P23 = (P2p + P3p) / 2;  P41 = (P4p + P1p) / 2
    i23_41 = (P23 - P41) / norm(P23 - P41)
    j1 = i23_41
    j2 = cross(j3, j1)

    T = hcat(j1, j2, j3)   # 3×3, columns are local axes
    return P0, T
end




"""
#########################################

MODEL generation

#########################################
"""



function collect_anchor(PLATE_table, SUPPORT_table, CONDSUP_table, tol)
    npl = length(PLATE_table)
    p_anchor = [zeros(0, 3) for _ in 1:npl]

    # anchor points from plate geometry intersections
    for ipl in 1:npl
        Q = zeros(4, 3)
        Q[1, :] = PLATE_table[ipl][1:3]
        Q[2, :] = PLATE_table[ipl][4:6]
        Q[3, :] = PLATE_table[ipl][7:9]
        Q[4, :] = PLATE_table[ipl][10:12]
        for iside in 1:4
            i1 = iside; i2 = iside == 4 ? 1 : iside + 1
            Q1 = Q[i1, :]; Q2 = Q[i2, :]
            for ipl2 in 1:npl
                ipl2 == ipl && continue
                for icorner in 1:4
                    idx = (icorner-1)*3+1 : (icorner-1)*3+3
                    P = PLATE_table[ipl2][idx]
                    if abs(Q1[2]-Q2[2]) < tol && abs(Q1[3]-Q2[3]) < tol && P[1] > min(Q1[1],Q2[1])+tol && P[1] < max(Q1[1],Q2[1])-tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], [P[1] Q1[2] Q1[3]])
                    end
                    if abs(Q1[1]-Q2[1]) < tol && abs(Q1[3]-Q2[3]) < tol && P[2] > min(Q1[2],Q2[2])+tol && P[2] < max(Q1[2],Q2[2])-tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], [Q1[1] P[2] Q1[3]])
                    end
                    if abs(Q1[1]-Q2[1]) < tol && abs(Q1[2]-Q2[2]) < tol && P[3] > min(Q1[3],Q2[3])+tol && P[3] < max(Q1[3],Q2[3])-tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], [Q1[1] Q1[2] P[3]])
                    end
                end
            end
        end
    end

    # anchor points from ordinary supports
    for ipl in 1:npl
        Q = zeros(4, 3)
        Q[1, :] = PLATE_table[ipl][1:3]
        Q[2, :] = PLATE_table[ipl][4:6]
        Q[3, :] = PLATE_table[ipl][7:9]
        Q[4, :] = PLATE_table[ipl][10:12]
        for iside in 1:4
            i1 = iside; i2 = iside == 4 ? 1 : iside + 1
            Q1 = Q[i1, :]; Q2 = Q[i2, :]
            for isup in 1:length(SUPPORT_table)
                stype = SUPPORT_table[isup][1]
                if stype in ("point", "p-line", "p-surf", "line", "surf")
                    jpl = SUPPORT_table[isup][2]
                    P1j = PLATE_table[jpl][1:3]; P2j = PLATE_table[jpl][4:6]
                    P3j = PLATE_table[jpl][7:9]; P4j = PLATE_table[jpl][10:12]
                    _, T = plate_transform(P1j, P2j, P3j, P4j)
                    if stype in ("point",)
                        npoi = 1; indpoi = [(3,4)]
                    else
                        npoi = 2; indpoi = [(3,4), (5,6)]
                    end
                    for ipoi in 1:npoi
                        ii1, ii2 = indpoi[ipoi]
                        Pp = [Float64(SUPPORT_table[isup][ii1]), Float64(SUPPORT_table[isup][ii2]), 0.0]
                        P = T * Pp + P1j
                        if abs(Q1[2]-Q2[2]) < tol && abs(Q1[3]-Q2[3]) < tol && P[1] > min(Q1[1],Q2[1])+tol && P[1] < max(Q1[1],Q2[1])-tol
                            p_anchor[ipl] = vcat(p_anchor[ipl], [P[1] Q1[2] Q1[3]])
                        end
                        if abs(Q1[1]-Q2[1]) < tol && abs(Q1[3]-Q2[3]) < tol && P[2] > min(Q1[2],Q2[2])+tol && P[2] < max(Q1[2],Q2[2])-tol
                            p_anchor[ipl] = vcat(p_anchor[ipl], [Q1[1] P[2] Q1[3]])
                        end
                        if abs(Q1[1]-Q2[1]) < tol && abs(Q1[2]-Q2[2]) < tol && P[3] > min(Q1[3],Q2[3])+tol && P[3] < max(Q1[3],Q2[3])-tol
                            p_anchor[ipl] = vcat(p_anchor[ipl], [Q1[1] Q1[2] P[3]])
                        end
                    end
                end
            end
        end
    end

    # anchor points from conditional supports
    for ipl in 1:npl
        Q = zeros(4, 3)
        Q[1, :] = PLATE_table[ipl][1:3]
        Q[2, :] = PLATE_table[ipl][4:6]
        Q[3, :] = PLATE_table[ipl][7:9]
        Q[4, :] = PLATE_table[ipl][10:12]
        for iside in 1:4
            i1 = iside; i2 = iside == 4 ? 1 : iside + 1
            Q1 = Q[i1, :]; Q2 = Q[i2, :]
            for isup in 1:length(CONDSUP_table)
                stype = CONDSUP_table[isup][1]
                jpl = CONDSUP_table[isup][2]
                P1j = PLATE_table[jpl][1:3]; P2j = PLATE_table[jpl][4:6]
                P3j = PLATE_table[jpl][7:9]; P4j = PLATE_table[jpl][10:12]
                _, T = plate_transform(P1j, P2j, P3j, P4j)
                pts_to_check = Vector{Vector{Float64}}()
                if stype == "point"
                    push!(pts_to_check, [Float64(CONDSUP_table[isup][3]), Float64(CONDSUP_table[isup][4]), 0.0])
                elseif stype == "surf"
                    push!(pts_to_check, [Float64(CONDSUP_table[isup][3]), Float64(CONDSUP_table[isup][4]), 0.0])
                    push!(pts_to_check, [Float64(CONDSUP_table[isup][5]), Float64(CONDSUP_table[isup][6]), 0.0])
                end
                for Pp in pts_to_check
                    P = T * Pp + P1j
                    if abs(Q1[2]-Q2[2]) < tol && abs(Q1[3]-Q2[3]) < tol && P[1] > min(Q1[1],Q2[1])+tol && P[1] < max(Q1[1],Q2[1])-tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], [P[1] Q1[2] Q1[3]])
                    end
                    if abs(Q1[1]-Q2[1]) < tol && abs(Q1[3]-Q2[3]) < tol && P[2] > min(Q1[2],Q2[2])+tol && P[2] < max(Q1[2],Q2[2])-tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], [Q1[1] P[2] Q1[3]])
                    end
                    if abs(Q1[1]-Q2[1]) < tol && abs(Q1[2]-Q2[2]) < tol && P[3] > min(Q1[3],Q2[3])+tol && P[3] < max(Q1[3],Q2[3])-tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], [Q1[1] Q1[2] P[3]])
                    end
                end
            end
        end
    end

    return p_anchor
end



function mesh4_rect(P1, P2, P3, P4, p_anchor_mat, size_x, size_y, tol)
    nmin = 1
    P0, T = plate_transform(P1, P2, P3, P4)

    P1s = T' * (P1 - P0)
    P2s = T' * (P2 - P0)
    P3s = T' * (P3 - P0)
    P4s = T' * (P4 - P0)

    x_anch = [P1s[1], P3s[1]]
    y_anch = [P1s[2], P3s[2]]
    for i in 1:size(p_anchor_mat, 1)
        Ps = T' * (p_anchor_mat[i, :] - P0)
        if minimum(abs.(x_anch .- Ps[1])) > tol
            push!(x_anch, Ps[1])
        end
        if minimum(abs.(y_anch .- Ps[2])) > tol
            push!(y_anch, Ps[2])
        end
    end
    x_anch = sort(x_anch)
    y_anch = sort(y_anch)

    x_seed = [x_anch[1]]
    for i in 2:length(x_anch)
        nx = max(round(Int, (x_anch[i] - x_anch[i-1]) / size_x), nmin)
        dx = (x_anch[i] - x_anch[i-1]) / nx
        append!(x_seed, x_seed[end] .+ (1:nx) .* dx)
    end
    nrx = length(x_seed)

    y_seed = [y_anch[1]]
    for i in 2:length(y_anch)
        ny = max(round(Int, (y_anch[i] - y_anch[i-1]) / size_y), nmin)
        dy = (y_anch[i] - y_anch[i-1]) / ny
        append!(y_seed, y_seed[end] .+ (1:ny) .* dy)
    end
    nry = length(y_seed)

    a = x_seed[nrx] - x_seed[1]
    b = y_seed[nry] - y_seed[1]
    x0 = x_seed[1]; y0 = y_seed[1]

    nodes = zeros(nrx * nry, 3)
    nr = 0
    for j in 1:nry, i in 1:nrx
        x = (x_seed[i] - x0) / a
        y = (y_seed[j] - y0) / b
        w1 = (1-x)*(1-y); w2 = x*(1-y); w3 = x*y; w4 = (1-x)*y
        Ps = w1*P1s[1:3] + w2*P2s[1:3] + w3*P3s[1:3] + w4*P4s[1:3]
        nr += 1
        nodes[nr, :] = T * Ps + P0
    end

    nelx = nrx - 1; nely = nry - 1
    elems = zeros(Int, nelx * nely, 4)
    nr = 0
    for j in 1:nely, i in 1:nelx
        nr += 1
        elems[nr, 1] = nrx*(j-1) + i
        elems[nr, 2] = nrx*(j-1) + i + 1
        elems[nr, 3] = nrx*j     + i + 1
        elems[nr, 4] = nrx*j     + i
    end
    return nodes, elems
end



function HOLE_TABLE_entry(HOLE_table, ipl)
    if ipl > length(HOLE_table) || isnothing(HOLE_table[ipl])
        return nothing
    end
    entry = HOLE_table[ipl]
    if isnothing(entry) || (isa(entry, HoleData) && isempty(entry.holes))
        return nothing
    end
    return entry
end


function elem_divide_by_plate(pl_elems, pl_nodes, PLATE_table, HOLE_table, tol)
    for ipl in 1:length(pl_elems)
        HOLE_TABLE_ipl = HOLE_TABLE_entry(HOLE_table, ipl)
        isnothing(HOLE_TABLE_ipl) && continue

        C1 = PLATE_table[ipl][1:3]; C2 = PLATE_table[ipl][4:6]
        C3 = PLATE_table[ipl][7:9]; C4 = PLATE_table[ipl][10:12]
        _, T = plate_transform(C1, C2, C3, C4)

        el_nodes = pl_elems[ipl]
        nelems = size(el_nodes, 1)
        eln2 = zeros(Int, 0, 9)
        elim_nodes = Int[]

        for ihole in 1:length(HOLE_TABLE_ipl.holes)
            npoly = size(HOLE_TABLE_ipl.holes[ihole], 1) - 1
            for iel in 1:nelems
                found = false; ipoly = 0
                while !found && ipoly < npoly
                    ipoly += 1
                    nsec = 0
                    P1h = HOLE_TABLE_ipl.holes[ihole][ipoly,   1:2]
                    P2h = HOLE_TABLE_ipl.holes[ihole][ipoly+1, 1:2]
                    el_sec_points = zeros(4, 4)
                    for iside in 1:4
                        i1 = iside; i2 = iside == 4 ? 1 : iside + 1
                        Q1g = pl_nodes[ipl][el_nodes[iel, i1], 1:3]
                        Q2g = pl_nodes[ipl][el_nodes[iel, i2], 1:3]
                        Q1 = (T' * (Q1g - C1))[1:2]
                        Q2 = (T' * (Q2g - C1))[1:2]
                        A = [P2h[1]-P1h[1]  Q1[1]-Q2[1];
                             P2h[2]-P1h[2]  Q1[2]-Q2[2]]
                        b_vec = [Q1[1]-P1h[1]; Q1[2]-P1h[2]]
                        if abs(det(A)) > tol
                            x = A \ b_vec
                            if x[1] >= -tol && x[1] <= 1+tol && x[2] >= -tol && x[2] <= 1+tol
                                if (x[1] >= tol && x[1] <= 1-tol) || (x[2] >= tol && x[2] <= 1-tol)
                                    # mesh error - just note it
                                else
                                    nsec += 1
                                    el_sec_points[nsec, :] = [ipoly, iel, iside, round(x[2])]
                                end
                            end
                        end
                    end
                    if nsec == 4
                        new_row = zeros(Int, 1, 9)
                        if el_sec_points[1, 4] == 1
                            new_row[1, 1:3] = el_nodes[iel, [3, 4, 2]]
                            el_nodes[iel, 1:3] = el_nodes[iel, [1, 2, 4]]
                        else
                            new_row[1, 1:3] = el_nodes[iel, [4, 1, 3]]
                            el_nodes[iel, 1:3] = el_nodes[iel, [2, 3, 1]]
                        end
                        eln2 = vcat(eln2, new_row)
                        push!(elim_nodes, iel)
                        found = true
                    end
                end
            end
        end

        for iel in elim_nodes
            el_nodes[iel, 4] = 0
        end
        pl_elems[ipl] = vcat(el_nodes, eln2)
    end
    return pl_elems
end


function four2three_by_plate(pl_elems, rat_q2t)
    for ipl in 1:length(pl_elems)
        el_nodes = pl_elems[ipl]
        nelems = size(el_nodes, 1)
        eln2 = zeros(Int, 0, 9)
        for iel in 1:nelems
            if el_nodes[iel, 4] > 0 && rand() > rat_q2t[ipl]
                new_row = zeros(Int, 1, 9)
                if rand() > (1 - rat_q2t[ipl]) / 2 + rat_q2t[ipl]
                    new_row[1, 1:3] = el_nodes[iel, [3, 4, 2]]
                    el_nodes[iel, 1:3] = el_nodes[iel, [1, 2, 4]]
                else
                    new_row[1, 1:3] = el_nodes[iel, [4, 1, 3]]
                    el_nodes[iel, 1:3] = el_nodes[iel, [2, 3, 1]]
                end
                el_nodes[iel, 4] = 0
                eln2 = vcat(eln2, new_row)
            end
        end
        pl_elems[ipl] = vcat(el_nodes, eln2)
    end
    return pl_elems
end



function edges_by_plate(pl_elems)
    pl_el_edges = Vector{Matrix{Int}}()
    pl_edge_nodes = Vector{Matrix{Int}}()

    for ipl in 1:length(pl_elems)
        el_nodes = pl_elems[ipl]
        nelems = size(el_nodes, 1)
        edge_tot = 0
        edge_nodes = zeros(Int, 0, 2)
        el_edges = zeros(Int, nelems, 4)

        for iel in 1:nelems
            nedg = el_nodes[iel, 4] == 0 ? 3 : 4
            for iedg in 1:nedg
                node1 = el_nodes[iel, iedg]
                node2 = iedg < nedg ? el_nodes[iel, iedg+1] : el_nodes[iel, 1]
                ind = findall(i -> (edge_nodes[i,1]==node1 && edge_nodes[i,2]==node2) ||
                                   (edge_nodes[i,1]==node2 && edge_nodes[i,2]==node1),
                              1:size(edge_nodes,1))
                if !isempty(ind)
                    el_edges[iel, iedg] = ind[1]
                else
                    edge_tot += 1
                    el_edges[iel, iedg] = edge_tot
                    edge_nodes = vcat(edge_nodes, [node1 node2])
                end
            end
        end
        push!(pl_el_edges, el_edges)
        push!(pl_edge_nodes, edge_nodes)
    end
    return pl_el_edges, pl_edge_nodes
end



function main_modelgen(PLATE_table, LOAD_table, SUPPORT_table, CONDSUP_table, HOLE_table,
                       FE_dat, fe_size, mesh_typ, rat_q2t, tolers)
    println("Model generation ")
    t0 = time()
    tol = tolers.zero_dist

    p_anchor = collect_anchor(PLATE_table, SUPPORT_table, CONDSUP_table, tol)

    pl_nodes = Vector{Matrix{Float64}}()
    pl_elems = Vector{Matrix{Int}}()

    for ipl in 1:length(PLATE_table)
        P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
        P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
        anch = isempty(p_anchor[ipl]) ? zeros(0, 3) : p_anchor[ipl]
        nod, ele = mesh4_rect(P1, P2, P3, P4, anch, fe_size[1], fe_size[2], tol)
        push!(pl_nodes, nod)
        nel = size(ele, 1)
        el_mat = zeros(Int, nel, 9)
        el_mat[:, 1:4] = ele
        push!(pl_elems, el_mat)
    end

    pl_elems = elem_divide_by_plate(pl_elems, pl_nodes, PLATE_table, HOLE_table, tol)
    pl_elems = four2three_by_plate(pl_elems, rat_q2t)
    pl_el_edges, pl_edge_nodes = edges_by_plate(pl_elems)

    # add global offsets to edge nodes
    nod_offsets = cumsum([0; [size(pl_nodes[i], 1) for i in 1:length(pl_nodes)-1]])
    for ipl in 1:length(PLATE_table)
        nod0 = nod_offsets[ipl]
        rows_en, _ = size(pl_edge_nodes[ipl])
        glo_edge = pl_edge_nodes[ipl] .+ nod0
        pl_edge_nodes[ipl] = hcat(pl_edge_nodes[ipl], glo_edge)  # cols 3,4 = global
    end

    # count total nodes and elements
    total_ele = sum(size(pl_elems[i], 1) for i in 1:length(pl_elems))
    total_nod = sum(size(pl_nodes[i], 1) for i in 1:length(pl_nodes))

    el_nodes = zeros(Int, total_ele, 9)
    el_props = zeros(total_ele, 5)
    nodes    = zeros(total_nod, 5)

    ele0 = 0; nod0 = 0
    for ipl in 1:length(PLATE_table)
        ele = size(pl_elems[ipl], 1)
        nod = size(pl_nodes[ipl], 1)
        nn = 9
        for iel in 1:ele
            ind0 = findall(j -> pl_elems[ipl][iel, j] == 0, 1:nn)
            base = fill(nod0, nn)
            base[ind0] .= 0
            el_nodes[ele0+iel, :] = pl_elems[ipl][iel, 1:nn] .+ base
        end
        el_props[ele0+1:ele0+ele, 1] .= ipl

        # update global edge node numbers
        pl_edge_nodes[ipl][:, 3:4] = pl_edge_nodes[ipl][:, 1:2] .+ nod0

        # assign element types
        ind_tri   = findall(j -> pl_elems[ipl][j, 4] == 0, 1:ele)
        ind_tri1  = findall(j -> pl_elems[ipl][j, 6] > 0 && pl_elems[ipl][j, 7] == 0, 1:ele)
        ind_quad  = findall(j -> pl_elems[ipl][j, 4] > 0 && pl_elems[ipl][j, 5] == 0, 1:ele)
        ind_quad2 = findall(j -> pl_elems[ipl][j, 8] > 0 && pl_elems[ipl][j, 9] == 0, 1:ele)
        el_props[ind_tri  .+ ele0, 2] .= mesh_typ[ipl, 1]
        el_props[ind_tri1 .+ ele0, 2] .= mesh_typ[ipl, 1]
        el_props[ind_quad .+ ele0, 2] .= mesh_typ[ipl, 2]
        el_props[ind_quad2.+ ele0, 2] .= mesh_typ[ipl, 2]

        el_props[ele0+1:ele0+ele, 3] .= PLATE_table[ipl][13]
        el_props[ele0+1:ele0+ele, 4] .= PLATE_table[ipl][14]
        nodes[nod0+1:nod0+nod, 1:3] = pl_nodes[ipl]
        nodes[nod0+1:nod0+nod, 4]   = collect(nod0+1 : nod0+nod)
        ele0 += ele; nod0 += nod
    end

    # hole markers
    el_props[:, 5] .= 1
    for ipl in 1:length(PLATE_table)
        HOLE_TABLE_ipl = HOLE_TABLE_entry(HOLE_table, ipl)
        isnothing(HOLE_TABLE_ipl) && continue
        indpl = findall(j -> el_props[j, 1] == ipl, 1:size(el_props, 1))
        P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
        P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
        _, T = plate_transform(P1, P2, P3, P4)
        for ihole in 1:length(HOLE_TABLE_ipl.holes)
            x0 = HOLE_TABLE_ipl.coord[ihole][1]
            y0 = HOLE_TABLE_ipl.coord[ihole][2]
            R  = HOLE_TABLE_ipl.radius[ihole]
            sec = HOLE_TABLE_ipl.sector[ihole]
            dalp2 = pi / sec
            R_ = R * cos(dalp2)
            C = T * [x0, y0, 0.0] + P1
            for iel in 1:size(el_nodes, 1)
                fe_id = Int(el_props[iel, 2])
                fe_id == 0 && continue
                elemnodenr = FE_dat[fe_id].nodenr
                midpoint = vec(mean(nodes[el_nodes[iel, 1:elemnodenr], 1:3], dims=1))
                dist = norm(midpoint - C)
                if dist < R_ + tol
                    el_props[iel, 5] = -1
                    el_props[iel, 3] /= 100
                end
            end
        end
    end

    # find duplicate nodes
    for in1 in 1:size(nodes, 1)-1
        nodes[in1, 4] != in1 && continue
        cands = findall(j -> abs(nodes[j,1]-nodes[in1,1]) < tol &&
                             abs(nodes[j,2]-nodes[in1,2]) < tol &&
                             abs(nodes[j,3]-nodes[in1,3]) < tol,
                        in1+1:size(nodes,1))
        for cin2 in cands
            in2 = in1 + cin2
            nodes[in2, 4] = 0
            nodes[in2, 5] = in1
        end
    end

    nodenr = count(j -> nodes[j, 4] > 0, 1:size(nodes, 1))
    nr = 0
    for i in 1:size(nodes, 1)
        if nodes[i, 4] > 0
            nr += 1
            nodes[i, 5] = nr
        end
    end
    for i in 1:size(nodes, 1)
        if nodes[i, 4] == 0
            nod = round(Int, nodes[i, 5])
            nodes[i, 5] = nodes[nod, 5]
        end
    end

    println("  done in $(round(time()-t0, digits=2)) s")
    return nodes, Int.(el_nodes), el_props, nodenr, pl_el_edges, pl_edge_nodes
end





"""
#########################################

LOAD VECTOR

#########################################
"""


function lv_el3_edge(elems_row, nodes, px, py, pz, dpn)
    T, P1e, P2e, P3e = ct_3node_g2e(elems_row, nodes)
    xe1=P1e[1]; ye1=P1e[2]
    xe2=P2e[1]; ye2=P2e[2]
    xe3=P3e[1]; ye3=P3e[2]

    L12 = sqrt((xe1-xe2)^2 + (ye1-ye2)^2)
    L23 = sqrt((xe2-xe3)^2 + (ye2-ye3)^2)
    L13 = sqrt((xe1-xe3)^2 + (ye1-ye3)^2)

    lv0_all = [
        [L12/3, L12/6, 0.0],
        [L12/6, L12/3, 0.0],
        [0.0, L23/3, L23/6],
        [0.0, L23/6, L23/3],
        [L13/6, 0.0, L13/3],
        [L13/3, 0.0, L13/6],
    ]

    elemnode = 3
    lv = zeros(elemnode * dpn)
    p0 = zeros(dpn)
    for i in 1:2*elemnode
        lv0 = lv0_all[i]
        pe = zeros(dpn)
        pe[1:min(3,dpn)] = [px[i], py[i], pz[i]][1:min(3,dpn)]
        pp = zeros(dpn*elemnode, elemnode)
        for k in 1:elemnode
            pp[(k-1)*dpn+1:k*dpn, k] = pe
        end
        lv += pp * lv0
    end
    return lv
end




function lv_el4_edge(elems_row, nodes, px, py, pz, dpn)
    T, P1e, P2e, P3e, P4e = ct_4node_g2e(elems_row, nodes)
    xe1=P1e[1]; ye1=P1e[2]
    xe2=P2e[1]; ye2=P2e[2]
    xe3=P3e[1]; ye3=P3e[2]
    xe4=P4e[1]; ye4=P4e[2]

    L12 = sqrt((xe1-xe2)^2 + (ye1-ye2)^2)
    L23 = sqrt((xe2-xe3)^2 + (ye2-ye3)^2)
    L34 = sqrt((xe3-xe4)^2 + (ye3-ye4)^2)
    L14 = sqrt((xe1-xe4)^2 + (ye1-ye4)^2)

    lv0_all = [
        [L12/3, L12/6, 0.0, 0.0],
        [L12/6, L12/3, 0.0, 0.0],
        [0.0, L23/3, L23/6, 0.0],
        [0.0, L23/6, L23/3, 0.0],
        [0.0, 0.0, L34/3, L34/6],
        [0.0, 0.0, L34/6, L34/3],
        [L14/3, 0.0, 0.0, L14/6],
        [L14/6, 0.0, 0.0, L14/3],
    ]

    elemnode = 4
    lv = zeros(elemnode * dpn)
    for i in 1:2*elemnode
        lv0 = lv0_all[i]
        pe = zeros(dpn)
        pe[1:min(3,dpn)] = [px[i], py[i], pz[i]][1:min(3,dpn)]
        pp = zeros(dpn*elemnode, elemnode)
        for k in 1:elemnode
            pp[(k-1)*dpn+1:k*dpn, k] = pe
        end
        lv += pp * lv0
    end
    return lv
end



function load_vec_point(nodes, nodenr, PLATE_table, LOAD_table, dpn, tolers)
    tol = tolers.zero_dist
    totaldof = nodenr * dpn
    LV = zeros(totaldof)

    for i in 1:length(LOAD_table)
        LOAD_table[i][1] != "point" && continue
        ipl = LOAD_table[i][2]
        L1p = [Float64(LOAD_table[i][3]), Float64(LOAD_table[i][4]), 0.0]
        px = Float64(LOAD_table[i][5]); py = Float64(LOAD_table[i][6]); pz = Float64(LOAD_table[i][7])
        mx = Float64(LOAD_table[i][8]); my = Float64(LOAD_table[i][9]); mz = Float64(LOAD_table[i][10])
        P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
        P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
        _, T3 = plate_transform(P1, P2, P3, P4)
        L1 = T3 * L1p + P1
        pvec = [px, py, pz, mx, my, mz]
        ltype_flag = LOAD_table[i][11]
        if isa(ltype_flag, String) && startswith(ltype_flag, "s")
            T6 = rotate3d(T3, 1, 6)
            pvec = T6 * pvec
        end
        dif = nodes[:, 1:3] .- L1'
        ind = findall(j -> abs(dif[j,1]) < tol && abs(dif[j,2]) < tol && abs(dif[j,3]) < tol, 1:size(nodes,1))
        isempty(ind) && continue
        active_nod = nodes[ind[1], 5]
        LV[(active_nod-1)*dpn+1 : (active_nod-1)*dpn+dpn] .+= pvec[1:dpn]
    end
    return LV
end



function load_vec_UDL(nodes, el_nodes, el_props, nodenr, PLATE_table, LOAD_table, FE_dat, dpn)
    nonzer = sum(FE_dat[Int(el_props[ii, 2])].nonzer_lv for ii in 1:size(el_nodes, 1))
    rows = zeros(Int, nonzer); cols = zeros(Int, nonzer); vals = zeros(nonzer)
    i0 = 0
    p0 = zeros(dpn)

    for i in 1:length(LOAD_table)
        LOAD_table[i][1] != "UDL" && continue
        ipl = LOAD_table[i][2]
        px = Float64(LOAD_table[i][7]); py = Float64(LOAD_table[i][8]); pz = Float64(LOAD_table[i][9])
        pvec = [px, py, pz, 0.0, 0.0, 0.0]

        ltype_flag = LOAD_table[i][10]
        if isa(ltype_flag, String) && startswith(ltype_flag, "s")
            P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
            P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
            _, T3 = plate_transform(P1, P2, P3, P4)
            T6 = rotate3d(T3, 1, 6)
            pvec = T6 * pvec
        end
        p = pvec[1:dpn]

        pl_elem_ids = findall(j -> el_props[j, 1] == ipl, 1:size(el_props, 1))
        for iii in pl_elem_ids
            fe_id = Int(el_props[iii, 2])
            elemnodenr = FE_dat[fe_id].nodenr
            elems_row = el_nodes[iii, 1:elemnodenr]

            if elemnodenr == 3
                T, P1e, P2e, P3e = ct_3node_g2e(elems_row, nodes)
                xe1=P1e[1]; ye1=P1e[2]; xe2=P2e[1]; ye2=P2e[2]; xe3=P3e[1]; ye3=P3e[2]
                lv00 = (xe1*ye2 - xe2*ye1 - xe1*ye3 + xe3*ye1 + xe2*ye3 - xe3*ye2) / 6
                lv0 = fill(lv00, 3)
                pe = zeros(dpn); pe[1:3] = T' * p[1:3]
                pp = zeros(dpn*3, 3)
                for k in 1:3; pp[(k-1)*dpn+1:k*dpn, k] = pe; end
            else
                T, P1e, P2e, P3e, P4e = ct_4node_g2e(elems_row, nodes)
                xe1=P1e[1]; ye1=P1e[2]; xe2=P2e[1]; ye2=P2e[2]
                xe3=P3e[1]; ye3=P3e[2]; xe4=P4e[1]; ye4=P4e[2]
                lv0 = [
                    xe1*ye2/6-xe2*ye1/6-xe1*ye4/6+xe2*ye3/12-xe3*ye2/12+xe4*ye1/6+xe2*ye4/12-xe4*ye2/12+xe3*ye4/12-xe4*ye3/12,
                    xe1*ye2/6-xe2*ye1/6-xe1*ye3/12+xe3*ye1/12-xe1*ye4/12+xe2*ye3/6-xe3*ye2/6+xe4*ye1/12+xe3*ye4/12-xe4*ye3/12,
                    xe1*ye2/12-xe2*ye1/12-xe1*ye4/12+xe2*ye3/6-xe3*ye2/6+xe4*ye1/12-xe2*ye4/12+xe4*ye2/12+xe3*ye4/6-xe4*ye3/6,
                    xe1*ye2/12-xe2*ye1/12+xe1*ye3/12-xe3*ye1/12-xe1*ye4/6+xe2*ye3/12-xe3*ye2/12+xe4*ye1/6+xe3*ye4/6-xe4*ye3/6
                ]
                pe = zeros(dpn); pe[1:3] = T' * p[1:3]
                pp = zeros(dpn*4, 4)
                for k in 1:4; pp[(k-1)*dpn+1:k*dpn, k] = pe; end
            end

            lve = pp * lv0
            TT = rotate3d(T, elemnodenr, dpn)
            lve = TT * lve

            ind = Int[]
            for jj in 1:elemnodenr
                nod_a = Int(nodes[elems_row[jj], 5])
                append!(ind, (nod_a-1)*dpn+1 : nod_a*dpn)
            end
            for r in 1:length(lve)
                v = lve[r]
                abs(v) > 0.0 || continue
                i0 += 1
                rows[i0] = ind[r]; cols[i0] = 1; vals[i0] = v
            end
        end
    end

    totaldof = nodenr * dpn
    LV = zeros(totaldof)
    for k in 1:i0
        LV[rows[k]] += vals[k]
    end
    return LV
end



function load_vec_edge(nodes, el_nodes, el_props, pl_el_edges, pl_edge_nodes, nodenr, PLATE_table, LOAD_table, FE_dat, dpn, tolers)
    tol = tolers.zero_dist
    nonzer = sum(FE_dat[Int(el_props[ii, 2])].nonzer_lv for ii in 1:size(el_nodes, 1))
    rows = zeros(Int, nonzer); vals = zeros(nonzer)
    i0 = 0

    for i in 1:length(LOAD_table)
        LOAD_table[i][1] != "edge" && continue
        ipl = LOAD_table[i][2]
        L1p = [Float64(LOAD_table[i][3]), Float64(LOAD_table[i][4]), 0.0]
        L2p = [Float64(LOAD_table[i][5]), Float64(LOAD_table[i][6]), 0.0]
        pvec1 = [Float64(LOAD_table[i][7]), Float64(LOAD_table[i][8]), Float64(LOAD_table[i][9]), 0.0, 0.0, 0.0]
        pvec2 = [Float64(LOAD_table[i][10]), Float64(LOAD_table[i][11]), Float64(LOAD_table[i][12]), 0.0, 0.0, 0.0]

        P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
        P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
        _, T = plate_transform(P1, P2, P3, P4)

        ltype_flag = LOAD_table[i][13]
        if isa(ltype_flag, String) && startswith(ltype_flag, "s")
            T6 = rotate3d(T, 1, 6)
            pvec1 = T6 * pvec1
            pvec2 = T6 * pvec2
        end

        pl_elem_ids = findall(j -> el_props[j, 1] == ipl, 1:size(el_props, 1))

        for iii in 1:size(pl_edge_nodes[ipl], 1)
            Q1 = nodes[pl_edge_nodes[ipl][iii, 3], 1:3]
            Q2 = nodes[pl_edge_nodes[ipl][iii, 4], 1:3]
            Q1p = T' * (Q1 - P1); Q2p = T' * (Q2 - P1)

            ind_el = findall(r -> any(pl_el_edges[ipl][r, :] .== iii), 1:size(pl_el_edges[ipl], 1))
            isempty(ind_el) && continue
            iel = ind_el[1]
            iedge_arr = findall(c -> pl_el_edges[ipl][iel, c] == iii, 1:size(pl_el_edges[ipl], 2))
            isempty(iedge_arr) && continue
            iedge = iedge_arr[1]

            fe_id = el_props[pl_elem_ids[iel], 2]
            elemnodenr = FE_dat[fe_id].nodenr
            vpx = zeros(2*elemnodenr); vpy = zeros(2*elemnodenr); vpz = zeros(2*elemnodenr)

            on_line1 = abs((Q1p[1]-L2p[1])*(L1p[2]-L2p[2]) + (L1p[1]-L2p[1])*(L2p[2]-Q1p[2])) < tol
            on_line2 = abs((Q2p[1]-L2p[1])*(L1p[2]-L2p[2]) + (L1p[1]-L2p[1])*(L2p[2]-Q2p[2])) < tol
            if on_line1 && on_line2
                L12 = norm(L2p - L1p)
                t1 = L12 > 0 ? norm(Q1p - L1p) / L12 : 0.0
                t2 = L12 > 0 ? norm(Q2p - L1p) / L12 : 0.0
                p1 = pvec1 + (pvec2 - pvec1) * t1
                p2 = pvec1 + (pvec2 - pvec1) * t2
                vpx[2*iedge-1]=p1[1]; vpy[2*iedge-1]=p1[2]; vpz[2*iedge-1]=p1[3]
                vpx[2*iedge]=p2[1];   vpy[2*iedge]=p2[2];   vpz[2*iedge]=p2[3]
            end

            if (maximum(abs.(vpx)) + maximum(abs.(vpy)) + maximum(abs.(vpz))) > 1e-8
                act_row = pl_elem_ids[iel]
                elems_row = el_nodes[act_row, 1:elemnodenr]
                if elemnodenr == 3
                    lve = lv_el3_edge(elems_row, nodes, vpx, vpy, vpz, dpn)
                else
                    lve = lv_el4_edge(elems_row, nodes, vpx, vpy, vpz, dpn)
                end

                ind = Int[]
                for jj in 1:elemnodenr
                    nod_a = Int(nodes[elems_row[jj], 5])
                    append!(ind, (nod_a-1)*dpn+1 : nod_a*dpn)
                end
                for r in 1:length(lve)
                    abs(lve[r]) > 0.0 || continue
                    i0 += 1
                    rows[i0] = ind[r]; vals[i0] = lve[r]
                end
            end
        end
    end

    totaldof = nodenr * dpn
    #LV = zeros(totaldof)
    #for k in 1:i0
    #    LV[rows[k]] += vals[k]
    #end
    LV=sparsevec(rows[1:i0],vals[1:i0],totaldof)
    LVT=LV'
    LV=LVT'
    return LV
end








"""
#########################################

STIFFNESS MATRIX

#########################################
"""



function ke_uv_3n_nocondens_6dof(P1, P2, P3, t, Ex, Ey, nuxy, nuyx, G, mat_sym)
    xe1=P1[1]; ye1=P1[2]
    xe2=P2[1]; ye2=P2[2]
    xe3=P3[1]; ye3=P3[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)

    D = xe1*ye2 - xe2*ye1 - xe1*ye3 + xe3*ye1 + xe2*ye3 - xe3*ye2

    ke_uv = [
        (G*t*xe2^2 - 2*G*t*xe2*xe3 + G*t*xe3^2 + E11*t*ye2^2 - 2*E11*t*ye2*ye3 + E11*t*ye3^2)/(2*D)  -(t*(xe2-xe3)*(ye2-ye3)*(E12+G))/(2*D)  -(E11*t*ye3^2+G*t*xe3^2+G*t*xe1*xe2-G*t*xe1*xe3-G*t*xe2*xe3+E11*t*ye1*ye2-E11*t*ye1*ye3-E11*t*ye2*ye3)/(2*D)  (E12*t*xe1*ye2-E12*t*xe1*ye3-E12*t*xe3*ye2+E12*t*xe3*ye3+G*t*xe2*ye1-G*t*xe3*ye1-G*t*xe2*ye3+G*t*xe3*ye3)/(2*D)  -(E11*t*ye2^2+G*t*xe2^2-G*t*xe1*xe2+G*t*xe1*xe3-G*t*xe2*xe3-E11*t*ye1*ye2+E11*t*ye1*ye3-E11*t*ye2*ye3)/(2*D)  -(E12*t*xe1*ye2-E12*t*xe1*ye3-E12*t*xe2*ye2+E12*t*xe2*ye3+G*t*xe2*ye1-G*t*xe2*ye2-G*t*xe3*ye1+G*t*xe3*ye2)/(2*D)
        -(t*(xe2-xe3)*(ye2-ye3)*(E21+G))/(2*D)  (E22*t*xe2^2-2*E22*t*xe2*xe3+E22*t*xe3^2+G*t*ye2^2-2*G*t*ye2*ye3+G*t*ye3^2)/(2*D)  (E21*t*xe2*ye1-E21*t*xe3*ye1-E21*t*xe2*ye3+E21*t*xe3*ye3+G*t*xe1*ye2-G*t*xe1*ye3-G*t*xe3*ye2+G*t*xe3*ye3)/(2*D)  -(E22*t*xe3^2+G*t*ye3^2+E22*t*xe1*xe2-E22*t*xe1*xe3-E22*t*xe2*xe3+G*t*ye1*ye2-G*t*ye1*ye3-G*t*ye2*ye3)/(2*D)  -(E21*t*xe2*ye1-E21*t*xe2*ye2-E21*t*xe3*ye1+E21*t*xe3*ye2+G*t*xe1*ye2-G*t*xe1*ye3-G*t*xe2*ye2+G*t*xe2*ye3)/(2*D)  -(E22*t*xe2^2+G*t*ye2^2-E22*t*xe1*xe2+E22*t*xe1*xe3-E22*t*xe2*xe3-G*t*ye1*ye2+G*t*ye1*ye3-G*t*ye2*ye3)/(2*D)
        -(E11*t*ye3^2+G*t*xe3^2+G*t*xe1*xe2-G*t*xe1*xe3-G*t*xe2*xe3+E11*t*ye1*ye2-E11*t*ye1*ye3-E11*t*ye2*ye3)/(2*D)  (E12*t*xe2*ye1-E12*t*xe3*ye1-E12*t*xe2*ye3+E12*t*xe3*ye3+G*t*xe1*ye2-G*t*xe1*ye3-G*t*xe3*ye2+G*t*xe3*ye3)/(2*D)  (t*(G*xe1^2-2*G*xe1*xe3+G*xe3^2+E11*ye1^2-2*E11*ye1*ye3+E11*ye3^2))/(2*D)  -(t*(xe1-xe3)*(ye1-ye3)*(E12+G))/(2*D)  -(t*(E11*ye1^2+G*xe1^2-G*xe1*xe2-G*xe1*xe3+G*xe2*xe3-E11*ye1*ye2-E11*ye1*ye3+E11*ye2*ye3))/(2*D)  (t*(E12*xe1*ye1-E12*xe2*ye1-E12*xe1*ye3+E12*xe2*ye3+G*xe1*ye1-G*xe1*ye2-G*xe3*ye1+G*xe3*ye2))/(2*D)
        (E21*t*xe1*ye2-E21*t*xe1*ye3-E21*t*xe3*ye2+E21*t*xe3*ye3+G*t*xe2*ye1-G*t*xe3*ye1-G*t*xe2*ye3+G*t*xe3*ye3)/(2*D)  -(E22*t*xe3^2+G*t*ye3^2+E22*t*xe1*xe2-E22*t*xe1*xe3-E22*t*xe2*xe3+G*t*ye1*ye2-G*t*ye1*ye3-G*t*ye2*ye3)/(2*D)  -(t*(xe1-xe3)*(ye1-ye3)*(E21+G))/(2*D)  (t*(E22*xe1^2-2*E22*xe1*xe3+E22*xe3^2+G*ye1^2-2*G*ye1*ye3+G*ye3^2))/(2*D)  (t*(E21*xe1*ye1-E21*xe1*ye2-E21*xe3*ye1+E21*xe3*ye2+G*xe1*ye1-G*xe2*ye1-G*xe1*ye3+G*xe2*ye3))/(2*D)  -(t*(E22*xe1^2+G*ye1^2-E22*xe1*xe2-E22*xe1*xe3+E22*xe2*xe3-G*ye1*ye2-G*ye1*ye3+G*ye2*ye3))/(2*D)
        -(E11*t*ye2^2+G*t*xe2^2-G*t*xe1*xe2+G*t*xe1*xe3-G*t*xe2*xe3-E11*t*ye1*ye2+E11*t*ye1*ye3-E11*t*ye2*ye3)/(2*D)  -(E12*t*xe2*ye1-E12*t*xe2*ye2-E12*t*xe3*ye1+E12*t*xe3*ye2+G*t*xe1*ye2-G*t*xe1*ye3-G*t*xe2*ye2+G*t*xe2*ye3)/(2*D)  -(t*(E11*ye1^2+G*xe1^2-G*xe1*xe2-G*xe1*xe3+G*xe2*xe3-E11*ye1*ye2-E11*ye1*ye3+E11*ye2*ye3))/(2*D)  (t*(E12*xe1*ye1-E12*xe1*ye2-E12*xe3*ye1+E12*xe3*ye2+G*xe1*ye1-G*xe2*ye1-G*xe1*ye3+G*xe2*ye3))/(2*D)  (t*(G*xe1^2-2*G*xe1*xe2+G*xe2^2+E11*ye1^2-2*E11*ye1*ye2+E11*ye2^2))/(2*D)  -(t*(xe1-xe2)*(ye1-ye2)*(E12+G))/(2*D)
        -(E21*t*xe1*ye2-E21*t*xe1*ye3-E21*t*xe2*ye2+E21*t*xe2*ye3+G*t*xe2*ye1-G*t*xe2*ye2-G*t*xe3*ye1+G*t*xe3*ye2)/(2*D)  -(E22*t*xe2^2+G*t*ye2^2-E22*t*xe1*xe2+E22*t*xe1*xe3-E22*t*xe2*xe3-G*t*ye1*ye2+G*t*ye1*ye3-G*t*ye2*ye3)/(2*D)  (t*(E21*xe1*ye1-E21*xe2*ye1-E21*xe1*ye3+E21*xe2*ye3+G*xe1*ye1-G*xe1*ye2-G*xe3*ye1+G*xe3*ye2))/(2*D)  -(t*(E22*xe1^2+G*ye1^2-E22*xe1*xe2-E22*xe1*xe3+E22*xe2*xe3-G*ye1*ye2-G*ye1*ye3+G*ye2*ye3))/(2*D)  -(t*(xe1-xe2)*(ye1-ye2)*(E21+G))/(2*D)  (t*(E22*xe1^2-2*E22*xe1*xe2+E22*xe2^2+G*ye1^2-2*G*ye1*ye2+G*ye2^2))/(2*D)
    ]

    return ke_uv
end



function ke_uv_3n_condens_from_8to6dof(P1, P2, P3, t, Ex, Ey, nuxy, nuyx, G, mat_sym)
    xe1=P1[1]; ye1=P1[2]
    xe2=P2[1]; ye2=P2[2]
    xe3=P3[1]; ye3=P3[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)

    D = xe1*ye2 - xe2*ye1 - xe1*ye3 + xe3*ye1 + xe2*ye3 - xe3*ye2

    kuv = [
        (G*t*xe2^2-2*G*t*xe2*xe3+G*t*xe3^2+E11*t*ye2^2-2*E11*t*ye2*ye3+E11*t*ye3^2)/(2*D)  -(t*(xe2-xe3)*(ye2-ye3)*(E12+G))/(2*D)  -(E11*t*ye3^2+G*t*xe3^2+G*t*xe1*xe2-G*t*xe1*xe3-G*t*xe2*xe3+E11*t*ye1*ye2-E11*t*ye1*ye3-E11*t*ye2*ye3)/(2*D)  (E12*t*xe1*ye2-E12*t*xe1*ye3-E12*t*xe3*ye2+E12*t*xe3*ye3+G*t*xe2*ye1-G*t*xe3*ye1-G*t*xe2*ye3+G*t*xe3*ye3)/(2*D)  -(E11*t*ye2^2+G*t*xe2^2-G*t*xe1*xe2+G*t*xe1*xe3-G*t*xe2*xe3-E11*t*ye1*ye2+E11*t*ye1*ye3-E11*t*ye2*ye3)/(2*D)  -(E12*t*xe1*ye2-E12*t*xe1*ye3-E12*t*xe2*ye2+E12*t*xe2*ye3+G*t*xe2*ye1-G*t*xe2*ye2-G*t*xe3*ye1+G*t*xe3*ye2)/(2*D)  0.0  0.0
        -(t*(xe2-xe3)*(ye2-ye3)*(E21+G))/(2*D)  (E22*t*xe2^2-2*E22*t*xe2*xe3+E22*t*xe3^2+G*t*ye2^2-2*G*t*ye2*ye3+G*t*ye3^2)/(2*D)  (E21*t*xe2*ye1-E21*t*xe3*ye1-E21*t*xe2*ye3+E21*t*xe3*ye3+G*t*xe1*ye2-G*t*xe1*ye3-G*t*xe3*ye2+G*t*xe3*ye3)/(2*D)  -(E22*t*xe3^2+G*t*ye3^2+E22*t*xe1*xe2-E22*t*xe1*xe3-E22*t*xe2*xe3+G*t*ye1*ye2-G*t*ye1*ye3-G*t*ye2*ye3)/(2*D)  -(E21*t*xe2*ye1-E21*t*xe2*ye2-E21*t*xe3*ye1+E21*t*xe3*ye2+G*t*xe1*ye2-G*t*xe1*ye3-G*t*xe2*ye2+G*t*xe2*ye3)/(2*D)  -(E22*t*xe2^2+G*t*ye2^2-E22*t*xe1*xe2+E22*t*xe1*xe3-E22*t*xe2*xe3-G*t*ye1*ye2+G*t*ye1*ye3-G*t*ye2*ye3)/(2*D)  0.0  0.0
        -(E11*t*ye3^2+G*t*xe3^2+G*t*xe1*xe2-G*t*xe1*xe3-G*t*xe2*xe3+E11*t*ye1*ye2-E11*t*ye1*ye3-E11*t*ye2*ye3)/(2*D)  (E12*t*xe2*ye1-E12*t*xe3*ye1-E12*t*xe2*ye3+E12*t*xe3*ye3+G*t*xe1*ye2-G*t*xe1*ye3-G*t*xe3*ye2+G*t*xe3*ye3)/(2*D)  (t*(G*xe1^2-2*G*xe1*xe3+G*xe3^2+E11*ye1^2-2*E11*ye1*ye3+E11*ye3^2))/(2*D)  -(t*(xe1-xe3)*(ye1-ye3)*(E12+G))/(2*D)  -(t*(E11*ye1^2+G*xe1^2-G*xe1*xe2-G*xe1*xe3+G*xe2*xe3-E11*ye1*ye2-E11*ye1*ye3+E11*ye2*ye3))/(2*D)  (t*(E12*xe1*ye1-E12*xe2*ye1-E12*xe1*ye3+E12*xe2*ye3+G*xe1*ye1-G*xe1*ye2-G*xe3*ye1+G*xe3*ye2))/(2*D)  0.0  0.0
        (E21*t*xe1*ye2-E21*t*xe1*ye3-E21*t*xe3*ye2+E21*t*xe3*ye3+G*t*xe2*ye1-G*t*xe3*ye1-G*t*xe2*ye3+G*t*xe3*ye3)/(2*D)  -(E22*t*xe3^2+G*t*ye3^2+E22*t*xe1*xe2-E22*t*xe1*xe3-E22*t*xe2*xe3+G*t*ye1*ye2-G*t*ye1*ye3-G*t*ye2*ye3)/(2*D)  -(t*(xe1-xe3)*(ye1-ye3)*(E21+G))/(2*D)  (t*(E22*xe1^2-2*E22*xe1*xe3+E22*xe3^2+G*ye1^2-2*G*ye1*ye3+G*ye3^2))/(2*D)  (t*(E21*xe1*ye1-E21*xe1*ye2-E21*xe3*ye1+E21*xe3*ye2+G*xe1*ye1-G*xe2*ye1-G*xe1*ye3+G*xe2*ye3))/(2*D)  -(t*(E22*xe1^2+G*ye1^2-E22*xe1*xe2-E22*xe1*xe3+E22*xe2*xe3-G*ye1*ye2-G*ye1*ye3+G*ye2*ye3))/(2*D)  0.0  0.0
        -(E11*t*ye2^2+G*t*xe2^2-G*t*xe1*xe2+G*t*xe1*xe3-G*t*xe2*xe3-E11*t*ye1*ye2+E11*t*ye1*ye3-E11*t*ye2*ye3)/(2*D)  -(E12*t*xe2*ye1-E12*t*xe2*ye2-E12*t*xe3*ye1+E12*t*xe3*ye2+G*t*xe1*ye2-G*t*xe1*ye3-G*t*xe2*ye2+G*t*xe2*ye3)/(2*D)  -(t*(E11*ye1^2+G*xe1^2-G*xe1*xe2-G*xe1*xe3+G*xe2*xe3-E11*ye1*ye2-E11*ye1*ye3+E11*ye2*ye3))/(2*D)  (t*(E12*xe1*ye1-E12*xe1*ye2-E12*xe3*ye1+E12*xe3*ye2+G*xe1*ye1-G*xe2*ye1-G*xe1*ye3+G*xe2*ye3))/(2*D)  (t*(G*xe1^2-2*G*xe1*xe2+G*xe2^2+E11*ye1^2-2*E11*ye1*ye2+E11*ye2^2))/(2*D)  -(t*(xe1-xe2)*(ye1-ye2)*(E12+G))/(2*D)  0.0  0.0
        -(E21*t*xe1*ye2-E21*t*xe1*ye3-E21*t*xe2*ye2+E21*t*xe2*ye3+G*t*xe2*ye1-G*t*xe2*ye2-G*t*xe3*ye1+G*t*xe3*ye2)/(2*D)  -(E22*t*xe2^2+G*t*ye2^2-E22*t*xe1*xe2+E22*t*xe1*xe3-E22*t*xe2*xe3-G*t*ye1*ye2+G*t*ye1*ye3-G*t*ye2*ye3)/(2*D)  (t*(E21*xe1*ye1-E21*xe2*ye1-E21*xe1*ye3+E21*xe2*ye3+G*xe1*ye1-G*xe1*ye2-G*xe3*ye1+G*xe3*ye2))/(2*D)  -(t*(E22*xe1^2+G*ye1^2-E22*xe1*xe2-E22*xe1*xe3+E22*xe2*xe3-G*ye1*ye2-G*ye1*ye3+G*ye2*ye3))/(2*D)  -(t*(xe1-xe2)*(ye1-ye2)*(E21+G))/(2*D)  (t*(E22*xe1^2-2*E22*xe1*xe2+E22*xe2^2+G*ye1^2-2*G*ye1*ye2+G*ye2^2))/(2*D)  0.0  0.0
        0.0  0.0  0.0  0.0  0.0  0.0  (81*t*(G*xe1^2-G*xe1*xe2-G*xe1*xe3+G*xe2^2-G*xe2*xe3+G*xe3^2+E11*ye1^2-E11*ye1*ye2-E11*ye1*ye3+E11*ye2^2-E11*ye2*ye3+E11*ye3^2))/(20*D)  (81*t*(E12+G)*(xe1*ye2-2*xe1*ye1+xe2*ye1+xe1*ye3-2*xe2*ye2+xe3*ye1+xe2*ye3+xe3*ye2-2*xe3*ye3))/(40*D)
        0.0  0.0  0.0  0.0  0.0  0.0  (81*t*(E21+G)*(xe1*ye2-2*xe1*ye1+xe2*ye1+xe1*ye3-2*xe2*ye2+xe3*ye1+xe2*ye3+xe3*ye2-2*xe3*ye3))/(40*D)  (81*t*(E22*xe1^2-E22*xe1*xe2-E22*xe1*xe3+E22*xe2^2-E22*xe2*xe3+E22*xe3^2+G*ye1^2-G*ye1*ye2-G*ye1*ye3+G*ye2^2-G*ye2*ye3+G*ye3^2))/(20*D)
    ]

    # static condensation
    inda = 1:6; indi = 7:8
    ke_uv = kuv[inda, inda] - kuv[inda, indi] * (kuv[indi, indi] \ kuv[indi, inda])
    return ke_uv
end



function ke_wt_3n_condens_from_12to9dof_MIN(P1, P2, P3, t, Ex, Ey, nuxy, nuyx, G, mat_sym)
    xe1=P1[1]; ye1=P1[2]
    xe2=P2[1]; ye2=P2[2]
    xe3=P3[1]; ye3=P3[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)

    D = xe1*ye2 - xe2*ye1 - xe1*ye3 + xe3*ye1 + xe2*ye3 - xe3*ye2

    # shear part (12×12)
    kwts = [
        (5*G*t*xe2^2-10*G*t*xe2*xe3+5*G*t*xe3^2+5*G*t*ye2^2-10*G*t*ye2*ye3+5*G*t*ye3^2)/(12*D)  (5*G*t*(xe2-xe3))/36  (5*G*t*(ye2-ye3))/36  -(5*G*t*xe3^2+5*G*t*ye3^2+5*G*t*xe1*xe2-5*G*t*xe1*xe3-5*G*t*xe2*xe3+5*G*t*ye1*ye2-5*G*t*ye1*ye3-5*G*t*ye2*ye3)/(12*D)  (5*G*t*(xe2-xe3))/36  (5*G*t*(ye2-ye3))/36  -(5*G*t*xe2^2+5*G*t*ye2^2-5*G*t*xe1*xe2+5*G*t*xe1*xe3-5*G*t*xe2*xe3-5*G*t*ye1*ye2+5*G*t*ye1*ye3-5*G*t*ye2*ye3)/(12*D)  (5*G*t*(xe2-xe3))/36  (5*G*t*(ye2-ye3))/36  -(5*G*t*(xe1*xe2-xe1*xe3+xe2*xe3+ye1*ye2-ye1*ye3+ye2*ye3-xe2^2-ye2^2))/(9*D)  -(5*G*t*(xe2^2-2*xe2*xe3+xe3^2+ye2^2-2*ye2*ye3+ye3^2))/(9*D)  (5*G*t*(xe1*xe2-xe1*xe3-xe2*xe3+ye1*ye2-ye1*ye3-ye2*ye3+xe3^2+ye3^2))/(9*D)
        (5*G*t*(xe2-xe3))/36  (5*G*t*D)/72  0.0  -(5*G*t*(xe1-xe3))/36  (5*G*t*D)/144  0.0  (5*G*t*(xe1-xe2))/36  (5*G*t*D)/144  0.0  (5*G*t*(xe2-2*xe1+xe3))/36  -(5*G*t*(xe2-xe3))/36  -(5*G*t*(xe2-2*xe1+xe3))/36
        (5*G*t*(ye2-ye3))/36  0.0  (5*G*t*D)/72  -(5*G*t*(ye1-ye3))/36  0.0  (5*G*t*D)/144  (5*G*t*(ye1-ye2))/36  0.0  (5*G*t*D)/144  (5*G*t*(ye2-2*ye1+ye3))/36  -(5*G*t*(ye2-ye3))/36  -(5*G*t*(ye2-2*ye1+ye3))/36
        -(5*G*t*xe3^2+5*G*t*ye3^2+5*G*t*xe1*xe2-5*G*t*xe1*xe3-5*G*t*xe2*xe3+5*G*t*ye1*ye2-5*G*t*ye1*ye3-5*G*t*ye2*ye3)/(12*D)  -(5*G*t*(xe1-xe3))/36  -(5*G*t*(ye1-ye3))/36  (5*G*t*(xe1^2-2*xe1*xe3+xe3^2+ye1^2-2*ye1*ye3+ye3^2))/(12*D)  -(5*G*t*(xe1-xe3))/36  -(5*G*t*(ye1-ye3))/36  (5*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(12*D)  -(5*G*t*(xe1-xe3))/36  -(5*G*t*(ye1-ye3))/36  -(5*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(9*D)  (5*G*t*(xe1*xe2-xe1*xe3-xe2*xe3+ye1*ye2-ye1*ye3-ye2*ye3+xe3^2+ye3^2))/(9*D)  -(5*G*t*(xe1^2-2*xe1*xe3+xe3^2+ye1^2-2*ye1*ye3+ye3^2))/(9*D)
        (5*G*t*(xe2-xe3))/36  (5*G*t*D)/144  0.0  -(5*G*t*(xe1-xe3))/36  (5*G*t*D)/72  0.0  (5*G*t*(xe1-xe2))/36  (5*G*t*D)/144  0.0  -(5*G*t*(xe1-2*xe2+xe3))/36  (5*G*t*(xe1-2*xe2+xe3))/36  (5*G*t*(xe1-xe3))/36
        (5*G*t*(ye2-ye3))/36  0.0  (5*G*t*D)/144  -(5*G*t*(ye1-ye3))/36  0.0  (5*G*t*D)/72  (5*G*t*(ye1-ye2))/36  0.0  (5*G*t*D)/144  -(5*G*t*(ye1-2*ye2+ye3))/36  (5*G*t*(ye1-2*ye2+ye3))/36  (5*G*t*(ye1-ye3))/36
        -(5*G*t*xe2^2+5*G*t*ye2^2-5*G*t*xe1*xe2+5*G*t*xe1*xe3-5*G*t*xe2*xe3-5*G*t*ye1*ye2+5*G*t*ye1*ye3-5*G*t*ye2*ye3)/(12*D)  (5*G*t*(xe1-xe2))/36  (5*G*t*(ye1-ye2))/36  (5*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(12*D)  (5*G*t*(xe1-xe2))/36  (5*G*t*(ye1-ye2))/36  (5*G*t*(xe1^2-2*xe1*xe2+xe2^2+ye1^2-2*ye1*ye2+ye2^2))/(12*D)  (5*G*t*(xe1-xe2))/36  (5*G*t*(ye1-ye2))/36  -(5*G*t*(xe1^2-2*xe1*xe2+xe2^2+ye1^2-2*ye1*ye2+ye2^2))/(9*D)  -(5*G*t*(xe1*xe2-xe1*xe3+xe2*xe3+ye1*ye2-ye1*ye3+ye2*ye3-xe2^2-ye2^2))/(9*D)  -(5*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(9*D)
        (5*G*t*(xe2-xe3))/36  (5*G*t*D)/144  0.0  -(5*G*t*(xe1-xe3))/36  (5*G*t*D)/144  0.0  (5*G*t*(xe1-xe2))/36  (5*G*t*D)/72  0.0  -(5*G*t*(xe1-xe2))/36  -(5*G*t*(xe1+xe2-2*xe3))/36  (5*G*t*(xe1+xe2-2*xe3))/36
        (5*G*t*(ye2-ye3))/36  0.0  (5*G*t*D)/144  -(5*G*t*(ye1-ye3))/36  0.0  (5*G*t*D)/144  (5*G*t*(ye1-ye2))/36  0.0  (5*G*t*D)/72  -(5*G*t*(ye1-ye2))/36  -(5*G*t*(ye1+ye2-2*ye3))/36  (5*G*t*(ye1+ye2-2*ye3))/36
        -(5*G*t*(xe1*xe2-xe1*xe3+xe2*xe3+ye1*ye2-ye1*ye3+ye2*ye3-xe2^2-ye2^2))/(9*D)  (5*G*t*(xe2-2*xe1+xe3))/36  (5*G*t*(ye2-2*ye1+ye3))/36  -(5*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(9*D)  -(5*G*t*(xe1-2*xe2+xe3))/36  -(5*G*t*(ye1-2*ye2+ye3))/36  -(5*G*t*(xe1^2-2*xe1*xe2+xe2^2+ye1^2-2*ye1*ye2+ye2^2))/(9*D)  -(5*G*t*(xe1-xe2))/36  -(5*G*t*(ye1-ye2))/36  -(10*G*t*(-xe1^2+xe1*xe2+xe1*xe3-xe2^2+xe2*xe3-xe3^2-ye1^2+ye1*ye2+ye1*ye3-ye2^2+ye2*ye3-ye3^2))/(9*D)  (10*G*t*(xe1*xe2-xe1*xe3+xe2*xe3+ye1*ye2-ye1*ye3+ye2*ye3-xe2^2-ye2^2))/(9*D)  (10*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(9*D)
        -(5*G*t*(xe2^2-2*xe2*xe3+xe3^2+ye2^2-2*ye2*ye3+ye3^2))/(9*D)  -(5*G*t*(xe2-xe3))/36  -(5*G*t*(ye2-ye3))/36  (5*G*t*(xe1*xe2-xe1*xe3-xe2*xe3+ye1*ye2-ye1*ye3-ye2*ye3+xe3^2+ye3^2))/(9*D)  (5*G*t*(xe1-2*xe2+xe3))/36  (5*G*t*(ye1-2*ye2+ye3))/36  -(5*G*t*(xe1*xe2-xe1*xe3+xe2*xe3+ye1*ye2-ye1*ye3+ye2*ye3-xe2^2-ye2^2))/(9*D)  -(5*G*t*(xe1+xe2-2*xe3))/36  -(5*G*t*(ye1+ye2-2*ye3))/36  (10*G*t*(xe1*xe2-xe1*xe3+xe2*xe3+ye1*ye2-ye1*ye3+ye2*ye3-xe2^2-ye2^2))/(9*D)  -(10*G*t*(-xe1^2+xe1*xe2+xe1*xe3-xe2^2+xe2*xe3-xe3^2-ye1^2+ye1*ye2+ye1*ye3-ye2^2+ye2*ye3-ye3^2))/(9*D)  -(10*G*t*(xe1*xe2-xe1*xe3-xe2*xe3+ye1*ye2-ye1*ye3-ye2*ye3+xe3^2+ye3^2))/(9*D)
        (5*G*t*(xe1*xe2-xe1*xe3-xe2*xe3+ye1*ye2-ye1*ye3-ye2*ye3+xe3^2+ye3^2))/(9*D)  -(5*G*t*(xe2-2*xe1+xe3))/36  -(5*G*t*(ye2-2*ye1+ye3))/36  -(5*G*t*(xe1^2-2*xe1*xe3+xe3^2+ye1^2-2*ye1*ye3+ye3^2))/(9*D)  (5*G*t*(xe1-xe3))/36  (5*G*t*(ye1-ye3))/36  -(5*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(9*D)  (5*G*t*(xe1+xe2-2*xe3))/36  (5*G*t*(ye1+ye2-2*ye3))/36  (10*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(9*D)  -(10*G*t*(xe1*xe2-xe1*xe3-xe2*xe3+ye1*ye2-ye1*ye3-ye2*ye3+xe3^2+ye3^2))/(9*D)  -(10*G*t*(-xe1^2+xe1*xe2+xe1*xe3-xe2^2+xe2*xe3-xe3^2-ye1^2+ye1*ye2+ye1*ye3-ye2^2+ye2*ye3-ye3^2))/(9*D)
    ]

    # bending part (12×12) - only non-zero entries at rotation DOFs (rows/cols 2,3,5,6,8,9)
    kwtb = zeros(12, 12)
    # row 2 entries
    kwtb[2,2]=(t^3*(E22*xe2^2-2*E22*xe2*xe3+E22*xe3^2+G*ye2^2-2*G*ye2*ye3+G*ye3^2))/(24*D)
    kwtb[2,3]=(t^3*(xe2-xe3)*(ye2-ye3)*(E21+G))/(24*D)
    kwtb[2,5]=-(t^3*(E22*xe3^2+G*ye3^2+E22*xe1*xe2-E22*xe1*xe3-E22*xe2*xe3+G*ye1*ye2-G*ye1*ye3-G*ye2*ye3))/(24*D)
    kwtb[2,6]=-(t^3*(E21*xe2*ye1-E21*xe3*ye1-E21*xe2*ye3+E21*xe3*ye3+G*xe1*ye2-G*xe1*ye3-G*xe3*ye2+G*xe3*ye3))/(24*D)
    kwtb[2,8]=-(t^3*(E22*xe2^2+G*ye2^2-E22*xe1*xe2+E22*xe1*xe3-E22*xe2*xe3-G*ye1*ye2+G*ye1*ye3-G*ye2*ye3))/(24*D)
    kwtb[2,9]= (t^3*(E21*xe2*ye1-E21*xe2*ye2-E21*xe3*ye1+E21*xe3*ye2+G*xe1*ye2-G*xe1*ye3-G*xe2*ye2+G*xe2*ye3))/(24*D)
    # row 3
    kwtb[3,2]=kwtb[2,3]
    kwtb[3,3]=(t^3*(G*xe2^2-2*G*xe2*xe3+G*xe3^2+E11*ye2^2-2*E11*ye2*ye3+E11*ye3^2))/(24*D)
    kwtb[3,5]=-(t^3*(E12*xe1*ye2-E12*xe1*ye3-E12*xe3*ye2+E12*xe3*ye3+G*xe2*ye1-G*xe3*ye1-G*xe2*ye3+G*xe3*ye3))/(24*D)
    kwtb[3,6]=-(t^3*(E11*ye3^2+G*xe3^2+G*xe1*xe2-G*xe1*xe3-G*xe2*xe3+E11*ye1*ye2-E11*ye1*ye3-E11*ye2*ye3))/(24*D)
    kwtb[3,8]= (t^3*(E12*xe1*ye2-E12*xe1*ye3-E12*xe2*ye2+E12*xe2*ye3+G*xe2*ye1-G*xe2*ye2-G*xe3*ye1+G*xe3*ye2))/(24*D)
    kwtb[3,9]=-(t^3*(E11*ye2^2+G*xe2^2-G*xe1*xe2+G*xe1*xe3-G*xe2*xe3-E11*ye1*ye2+E11*ye1*ye3-E11*ye2*ye3))/(24*D)
    # row 5
    kwtb[5,2]=kwtb[2,5]; kwtb[5,3]=kwtb[3,5]
    kwtb[5,5]=(t^3*(E22*xe1^2-2*E22*xe1*xe3+E22*xe3^2+G*ye1^2-2*G*ye1*ye3+G*ye3^2))/(24*D)
    kwtb[5,6]=(t^3*(xe1-xe3)*(ye1-ye3)*(E21+G))/(24*D)
    kwtb[5,8]=-(t^3*(E22*xe1^2+G*ye1^2-E22*xe1*xe2-E22*xe1*xe3+E22*xe2*xe3-G*ye1*ye2-G*ye1*ye3+G*ye2*ye3))/(24*D)
    kwtb[5,9]=-(t^3*(E21*xe1*ye1-E21*xe1*ye2-E21*xe3*ye1+E21*xe3*ye2+G*xe1*ye1-G*xe2*ye1-G*xe1*ye3+G*xe2*ye3))/(24*D)
    # row 6
    kwtb[6,2]=kwtb[2,6]; kwtb[6,3]=kwtb[3,6]; kwtb[6,5]=kwtb[5,6]
    kwtb[6,6]=(t^3*(G*xe1^2-2*G*xe1*xe3+G*xe3^2+E11*ye1^2-2*E11*ye1*ye3+E11*ye3^2))/(24*D)
    kwtb[6,8]=-(t^3*(E12*xe1*ye1-E12*xe2*ye1-E12*xe1*ye3+E12*xe2*ye3+G*xe1*ye1-G*xe1*ye2-G*xe3*ye1+G*xe3*ye2))/(24*D)
    kwtb[6,9]=-(t^3*(E11*ye1^2+G*xe1^2-G*xe1*xe2-G*xe1*xe3+G*xe2*xe3-E11*ye1*ye2-E11*ye1*ye3+E11*ye2*ye3))/(24*D)
    # row 8
    kwtb[8,2]=kwtb[2,8]; kwtb[8,3]=kwtb[3,8]; kwtb[8,5]=kwtb[5,8]; kwtb[8,6]=kwtb[6,8]
    kwtb[8,8]=(t^3*(E22*xe1^2-2*E22*xe1*xe2+E22*xe2^2+G*ye1^2-2*G*ye1*ye2+G*ye2^2))/(24*D)
    kwtb[8,9]=(t^3*(xe1-xe2)*(ye1-ye2)*(E21+G))/(24*D)
    # row 9
    kwtb[9,2]=kwtb[2,9]; kwtb[9,3]=kwtb[3,9]; kwtb[9,5]=kwtb[5,9]; kwtb[9,6]=kwtb[6,9]; kwtb[9,8]=kwtb[8,9]
    kwtb[9,9]=(t^3*(G*xe1^2-2*G*xe1*xe2+G*xe2^2+E11*ye1^2-2*E11*ye1*ye2+E11*ye2^2))/(24*D)

    # static condensation
    inda = 1:9; indi = 10:12
    kwts_c = kwts[inda, inda] - kwts[inda, indi] * (kwts[indi, indi] \ kwts[indi, inda])
    kwtb_c = kwtb[inda, inda]

    # shear correction
    alpha = sum(diag(kwts_c[4:9, 4:9])) / max(sum(diag(kwtb_c[4:9, 4:9])), 1e-30)
    Cs = 0.0
    ke_wt = (1 / (1 + Cs*alpha)) * kwts_c + kwtb_c
    return ke_wt
end



function ke_uv_4n_nocondens_8dof_num(P1, P2, P3, P4, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp_in)
    xe1=P1[1]; ye1=P1[2]
    xe2=P2[1]; ye2=P2[2]
    xe3=P3[1]; ye3=P3[2]
    xe4=P4[1]; ye4=P4[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)
    D = [E11 E12 0.0; E21 E22 0.0; 0.0 0.0 G]
    Duv = D * t

    nGp = round(Int, sqrt(nGp_in))
    locx, locy, wei = GL4(nGp)

    kuv = zeros(8, 8)
    for i in 1:nGp, j in 1:nGp
        x = locx[i, j]; y = locy[i, j]
        J = [xe1*(y/4-1/4)-xe2*(y/4-1/4)+xe3*(y/4+1/4)-xe4*(y/4+1/4)  ye1*(y/4-1/4)-ye2*(y/4-1/4)+ye3*(y/4+1/4)-ye4*(y/4+1/4);
             xe1*(x/4-1/4)-xe2*(x/4+1/4)+xe3*(x/4+1/4)-xe4*(x/4-1/4)  ye1*(x/4-1/4)-ye2*(x/4+1/4)+ye3*(x/4+1/4)-ye4*(x/4-1/4)]
        detJ = det(J); Ji = inv(J)
        Ji11=Ji[1,1]; Ji12=Ji[1,2]; Ji21=Ji[2,1]; Ji22=Ji[2,2]
        N1x=(y-1)/4; N2x=-N1x; N3x=(y+1)/4; N4x=-N3x
        N1y=(x-1)/4; N2y=-(x+1)/4; N3y=-N2y; N4y=-N1y
        Buv = [
            N1x*Ji11+N1y*Ji12  0.0  N2x*Ji11+N2y*Ji12  0.0  N3x*Ji11+N3y*Ji12  0.0  N4x*Ji11+N4y*Ji12  0.0;
            0.0  N1x*Ji21+N1y*Ji22  0.0  N2x*Ji21+N2y*Ji22  0.0  N3x*Ji21+N3y*Ji22  0.0  N4x*Ji21+N4y*Ji22;
            N1x*Ji21+N1y*Ji22  N1x*Ji11+N1y*Ji12  N2x*Ji21+N2y*Ji22  N2x*Ji11+N2y*Ji12  N3x*Ji21+N3y*Ji22  N3x*Ji11+N3y*Ji12  N4x*Ji21+N4y*Ji22  N4x*Ji11+N4y*Ji12
        ]
        kuv += Buv' * Duv * Buv * detJ * wei[i, j]
    end
    return kuv
end



function ke_uv_4n_condens_from_12to8dof_num(P1, P2, P3, P4, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp_in)
    xe1=P1[1]; ye1=P1[2]
    xe2=P2[1]; ye2=P2[2]
    xe3=P3[1]; ye3=P3[2]
    xe4=P4[1]; ye4=P4[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)
    D = [E11 E12 0.0; E21 E22 0.0; 0.0 0.0 G]
    Duv = D * t

    nGp = round(Int, sqrt(nGp_in))
    locx, locy, wei = GL4(nGp)

    kuv = zeros(12, 12)
    for i in 1:nGp, j in 1:nGp
        x = locx[i, j]; y = locy[i, j]
        J = [xe1*(y/4-1/4)-xe2*(y/4-1/4)+xe3*(y/4+1/4)-xe4*(y/4+1/4)  ye1*(y/4-1/4)-ye2*(y/4-1/4)+ye3*(y/4+1/4)-ye4*(y/4+1/4);
             xe1*(x/4-1/4)-xe2*(x/4+1/4)+xe3*(x/4+1/4)-xe4*(x/4-1/4)  ye1*(x/4-1/4)-ye2*(x/4+1/4)+ye3*(x/4+1/4)-ye4*(x/4-1/4)]
        detJ = det(J); Ji = inv(J)
        Ji11=Ji[1,1]; Ji12=Ji[1,2]; Ji21=Ji[2,1]; Ji22=Ji[2,2]
        N1x=(y-1)/4; N2x=-N1x; N3x=(y+1)/4; N4x=-N3x
        N5x=-2*x;   N6x=0.0
        N1y=(x-1)/4; N2y=-(x+1)/4; N3y=-N2y; N4y=-N1y
        N5y=0.0;    N6y=-2*y
        Buv = [
            N1x*Ji11+N1y*Ji12  0.0  N2x*Ji11+N2y*Ji12  0.0  N3x*Ji11+N3y*Ji12  0.0  N4x*Ji11+N4y*Ji12  0.0  N5x*Ji11+N5y*Ji12  0.0  N6x*Ji11+N6y*Ji12  0.0;
            0.0  N1x*Ji21+N1y*Ji22  0.0  N2x*Ji21+N2y*Ji22  0.0  N3x*Ji21+N3y*Ji22  0.0  N4x*Ji21+N4y*Ji22  0.0  N5x*Ji21+N5y*Ji22  0.0  N6x*Ji21+N6y*Ji22;
            N1x*Ji21+N1y*Ji22  N1x*Ji11+N1y*Ji12  N2x*Ji21+N2y*Ji22  N2x*Ji11+N2y*Ji12  N3x*Ji21+N3y*Ji22  N3x*Ji11+N3y*Ji12  N4x*Ji21+N4y*Ji22  N4x*Ji11+N4y*Ji12  N5x*Ji21+N5y*Ji22  N5x*Ji11+N5y*Ji12  N6x*Ji21+N6y*Ji22  N6x*Ji11+N6y*Ji12
        ]
        kuv += Buv' * Duv * Buv * detJ * wei[i, j]
    end

    inda = 1:8; indi = 9:12
    ke_uv = kuv[inda, inda] - kuv[inda, indi] * (kuv[indi, indi] \ kuv[indi, inda])
    return ke_uv
end



function ke_wt_4n_condens_from_18to12dof_num(P1, P2, P3, P4, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp_in)
    X1=P1[1]; Y1=P1[2]
    X2=P2[1]; Y2=P2[2]
    X3=P3[1]; Y3=P3[2]
    X4=P4[1]; Y4=P4[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)
    D   = [E11 E12 0.0; E21 E22 0.0; 0.0 0.0 G]
    Dwtb = D * t^3 / 12
    Dwts = Matrix{Float64}(I, 2, 2) * (5/6 * G * t)

    nGp = round(Int, sqrt(nGp_in))
    locx, locy, wei = GL4(nGp)

    kwts = zeros(18, 18)
    kwtb = zeros(18, 18)

    for i in 1:nGp, j in 1:nGp
        x = locx[i, j]; y = locy[i, j]
        J = [X1*(y/4-1/4)-X2*(y/4-1/4)+X3*(y/4+1/4)-X4*(y/4+1/4)  Y1*(y/4-1/4)-Y2*(y/4-1/4)+Y3*(y/4+1/4)-Y4*(y/4+1/4);
             X1*(x/4-1/4)-X2*(x/4+1/4)+X3*(x/4+1/4)-X4*(x/4-1/4)  Y1*(x/4-1/4)-Y2*(x/4+1/4)+Y3*(x/4+1/4)-Y4*(x/4-1/4)]
        detJ = det(J); Ji = inv(J)
        Ji11=Ji[1,1]; Ji12=Ji[1,2]; Ji21=Ji[2,1]; Ji22=Ji[2,2]
        N1=(x-1)*(y-1)/4; N2=-(x+1)*(y-1)/4; N3=(x+1)*(y+1)/4; N4=-(x-1)*(y+1)/4
        N5=1-x^2; N6=1-y^2
        N1x=(y-1)/4; N2x=-N1x; N3x=(y+1)/4; N4x=-N3x; N5x=-2*x; N6x=0.0
        N1y=(x-1)/4; N2y=-(x+1)/4; N3y=-N2y; N4y=-N1y; N5y=0.0; N6y=-2*y

        Bwts = [
            N1x*Ji11+N1y*Ji12  0.0  N1  N2x*Ji11+N2y*Ji12  0.0  N2  N3x*Ji11+N3y*Ji12  0.0  N3  N4x*Ji11+N4y*Ji12  0.0  N4  N5x*Ji11+N5y*Ji12  0.0  N5  N6x*Ji11+N6y*Ji12  0.0  N6;
            N1x*Ji21+N1y*Ji22  -N1  0.0  N2x*Ji21+N2y*Ji22  -N2  0.0  N3x*Ji21+N3y*Ji22  -N3  0.0  N4x*Ji21+N4y*Ji22  -N4  0.0  N5x*Ji21+N5y*Ji22  -N5  0.0  N6x*Ji21+N6y*Ji22  -N6  0.0
        ]
        kwts += Bwts' * Dwts * Bwts * detJ * wei[i, j]

        Bwtb = zeros(3, 18)
        # row 1: d(ty)/dx
        for k in 1:6
            Nkx = [N1x, N2x, N3x, N4x, N5x, N6x][k]
            Nky = [N1y, N2y, N3y, N4y, N5y, N6y][k]
            col = (k-1)*3 + 2
            Bwtb[1, col] = Nkx*Ji11 + Nky*Ji12
        end
        # row 2: -d(tx)/dy  (note tx is 3rd DOF per node)
        for k in 1:6
            Nkx = [N1x, N2x, N3x, N4x, N5x, N6x][k]
            Nky = [N1y, N2y, N3y, N4y, N5y, N6y][k]
            col = (k-1)*3 + 3
            Bwtb[2, col] = -(Nkx*Ji21 + Nky*Ji22)
        end
        # row 3: d(ty)/dy - d(tx)/dx  -> Bwtb[3,:] = -Bwtb[1,tx_cols] and Bwtb[3,ty_cols] refers to col+1
        for k in 1:6
            Nkx = [N1x, N2x, N3x, N4x, N5x, N6x][k]
            Nky = [N1y, N2y, N3y, N4y, N5y, N6y][k]
            col_ty = (k-1)*3 + 2
            col_tx = (k-1)*3 + 3
            Bwtb[3, col_ty] = Nkx*Ji21 + Nky*Ji22
            Bwtb[3, col_tx] = -(Nkx*Ji11 + Nky*Ji12)
        end
        kwtb += Bwtb' * Dwtb * Bwtb * detJ * wei[i, j]
    end

    inda = 1:12; indi = 13:18
    kwt1 = kwts + kwtb
    ke_wt = kwt1[inda, inda] - kwt1[inda, indi] * (kwt1[indi, indi] \ kwt1[indi, inda])
    return ke_wt
end



function add_drill(k5dpn::AbstractMatrix, nodenr::Int, ifelastic::Bool)
    ind = Int[]
    i0 = 0
    for i in 1:nodenr
        append!(ind, (1:5) .+ i0)
        i0 += 6
    end
    k6dpn = zeros(6 * nodenr, 6 * nodenr)
    k6dpn[ind, ind] = k5dpn

    if ifelastic
        ind2 = Int[]
        i0 = 0
        for i in 1:nodenr
            append!(ind2, [4, 5] .+ i0)
            i0 += 5
        end
        kd = diag(k5dpn)
        stif = minimum(kd[ind2]) / 100
        for i in 1:nodenr
            k6dpn[6*i, 6*i] = stif
        end
    end
    return k6dpn
end



function stiffmat_e(nodes, el_nodes, el_props, nodenr, MAT_table, FE_dat, mat_sym, dpn)
    nonzer = sum(FE_dat[Int(el_props[ii, 2])].nonzer_ke for ii in 1:size(el_nodes, 1))
    rows = zeros(Int, nonzer)
    cols = zeros(Int, nonzer)
    vals = zeros(nonzer)
    i0 = 0

    for ii in 1:size(el_nodes, 1)
        fe_id = Int(el_props[ii, 2])
        tel = el_props[ii, 3]
        imat = Int(el_props[ii, 4])
        mv = MAT_table[imat]
        Ex = mv[1]; Ey = mv[2]; nuxy = mv[3]; nuyx = mv[4]; G = mv[5]
        fd = FE_dat[fe_id]
        elemnodenr = fd.nodenr
        elems_row = el_nodes[ii, 1:elemnodenr]

        if elemnodenr == 3
            T, P1e, P2e, P3e = ct_3node_g2e(elems_row, nodes)
        else
            T, P1e, P2e, P3e, P4e = ct_4node_g2e(elems_row, nodes)
        end

        # membrane stiffness
        if elemnodenr == 3 && fd.m_dof_ini == 6 && !fd.m_ifcondens && !fd.m_numint
            ke_uv = ke_uv_3n_nocondens_6dof(P1e, P2e, P3e, tel, Ex, Ey, nuxy, nuyx, G, mat_sym)
        elseif elemnodenr == 3 && fd.m_dof_ini == 8 && fd.m_ifcondens && !fd.m_numint
            ke_uv = ke_uv_3n_condens_from_8to6dof(P1e, P2e, P3e, tel, Ex, Ey, nuxy, nuyx, G, mat_sym)
        elseif elemnodenr == 4 && fd.m_dof_ini == 8 && !fd.m_ifcondens && fd.m_numint
            ke_uv = ke_uv_4n_nocondens_8dof_num(P1e, P2e, P3e, P4e, tel, Ex, Ey, nuxy, nuyx, G, mat_sym, fd.m_nGpe)
        elseif elemnodenr == 4 && fd.m_dof_ini == 12 && fd.m_ifcondens && fd.m_numint
            ke_uv = ke_uv_4n_condens_from_12to8dof_num(P1e, P2e, P3e, P4e, tel, Ex, Ey, nuxy, nuyx, G, mat_sym, fd.m_nGpe)
        else
            ke_uv = zeros(elemnodenr*2, elemnodenr*2)
        end

        # bending stiffness
        if elemnodenr == 3 && fd.b_dof_ini == 12 && fd.b_ifcondens && !fd.b_numint
            ke_wt = ke_wt_3n_condens_from_12to9dof_MIN(P1e, P2e, P3e, tel, Ex, Ey, nuxy, nuyx, G, mat_sym)
        elseif elemnodenr == 4 && fd.b_dof_ini == 18 && fd.b_ifcondens && fd.b_numint
            ke_wt = ke_wt_4n_condens_from_18to12dof_num(P1e, P2e, P3e, P4e, tel, Ex, Ey, nuxy, nuyx, G, mat_sym, fd.b_nGpe)
        else
            ke_wt = zeros(elemnodenr*3, elemnodenr*3)
        end

        # assemble ke from membrane + bending parts (5 DOF per node)
        induv = Int[]; indwt = Int[]
        n0 = 0
        for i in 1:elemnodenr
            append!(induv, [n0+1, n0+2])
            append!(indwt, [n0+3, n0+4, n0+5])
            n0 += 5
        end
        ke = zeros(elemnodenr*5, elemnodenr*5)
        ke[induv, induv] = ke_uv
        ke[indwt, indwt] = ke_wt

        if dpn == 6
            ke = add_drill(ke, elemnodenr, true)
        end

        TT = rotate3d(T, elemnodenr, dpn)
        ke = TT * ke * TT'

        # find node global DOF indices
        ind = Int[]
        for jj in 1:elemnodenr
            nod = elems_row[jj]
            nod_active = Int(nodes[nod, 5])
            append!(ind, (nod_active-1)*dpn+1 : nod_active*dpn)
        end

        # store sparse entries
        for r in 1:size(ke, 1), c in 1:size(ke, 2)
            v = ke[r, c]
            if abs(v) > 0.0
                i0 += 1
                if i0 > length(rows)
                    resize!(rows, 2*i0); resize!(cols, 2*i0); resize!(vals, 2*i0)
                end
                rows[i0] = ind[r]
                cols[i0] = ind[c]
                vals[i0] = v
            end
        end
    end

    totaldof = nodenr * dpn
    Ke = sparse(rows[1:i0], cols[1:i0], vals[1:i0], totaldof, totaldof)
    return Ke
end




"""
#########################################

SUPPORTS

#########################################
"""


function fem_support(Ke, PLATE_table, SUPPORT_table, FE_dat, nodes, el_nodes, el_props, nodenr, dpn, tolers)
    tol = tolers.zero_dist
    zero_spring = tolers.zero_spring
    rigid = maximum(Ke) * 1000

    nonzer = nodenr * dpn^2
    rows = zeros(Int, nonzer); cols = zeros(Int, nonzer); vals = zeros(nonzer)
    i0 = 0
    nsup = length(SUPPORT_table)

    # point/line/surface-point supports
    for i in 1:nsup
        stype = SUPPORT_table[i][1]
        if stype in ("point", "p-line", "p-surf")
            ipl = SUPPORT_table[i][2]
            P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
            P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
            _, T = plate_transform(P1, P2, P3, P4)

            if stype == "point"
                Qp = [Float64(SUPPORT_table[i][3]), Float64(SUPPORT_table[i][4]), 0.0]
                Q = T * Qp + P1
                nods = findall(j -> abs(nodes[j,1]-Q[1]) < tol && abs(nodes[j,2]-Q[2]) < tol &&
                               abs(nodes[j,3]-Q[3]) < tol && nodes[j,4] > 0, 1:size(nodes,1))
                jadd = 0
            elseif stype == "p-line"
                Q1p = [Float64(SUPPORT_table[i][3]), Float64(SUPPORT_table[i][4]), 0.0]
                Q2p = [Float64(SUPPORT_table[i][5]), Float64(SUPPORT_table[i][6]), 0.0]
                nodesp = (nodes[:, 1:3] .- P1') * T
                nods = findall(j -> abs((nodesp[j,1]-Q2p[1])*(Q1p[2]-Q2p[2]) +
                               (Q1p[1]-Q2p[1])*(Q2p[2]-nodesp[j,2])) < tol &&
                               abs(nodesp[j,3]) < tol, 1:size(nodes,1))
                jadd = 2
            else  # p-surf
                Q1p = [Float64(SUPPORT_table[i][3]), Float64(SUPPORT_table[i][4]), 0.0]
                Q2p = [Float64(SUPPORT_table[i][5]), Float64(SUPPORT_table[i][6]), 0.0]
                nodesp = (nodes[:, 1:3] .- P1') * T
                nods = findall(j -> nodesp[j,1] >= Q1p[1]-tol && nodesp[j,2] >= Q1p[2]-tol &&
                               nodesp[j,3] >= Q1p[3]-tol && nodesp[j,1] <= Q2p[1]+tol &&
                               nodesp[j,2] <= Q2p[2]+tol && nodesp[j,3] <= Q2p[3]+tol &&
                               nodes[j,4] > tol, 1:size(nodes,1))
                jadd = 2
            end

            stif = zeros(6)
            for j in 1:6
                sv = Float64(SUPPORT_table[i][j+4+jadd])
                sv < -1e-6 && (stif[j] = rigid)
                sv > zero_spring && (stif[j] = sv)
            end
            stif_mat = diagm(0 => stif)
            coord_flag = SUPPORT_table[i][11+jadd]
            if isa(coord_flag, String) && startswith(coord_flag, "s")
                stif_mat[1:3, 1:3] = T * stif_mat[1:3, 1:3] * T'
                stif_mat[4:6, 4:6] = T * stif_mat[4:6, 4:6] * T'
            end

            for inod in nods
                nod_active = Int(nodes[inod, 5])
                ndof = (nod_active - 1) * dpn
                for idof in 1:dpn, jdof in 1:dpn
                    st = stif_mat[idof, jdof]
                    if abs(st) > zero_spring / 1000
                        i0 += 1
                        rows[i0] = ndof + idof; cols[i0] = ndof + jdof; vals[i0] = st
                    end
                end
            end
        end
    end

    # surface supports
    for i in 1:nsup
        SUPPORT_table[i][1] != "surf" && continue
        ipl = SUPPORT_table[i][2]
        P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
        P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
        _, T = plate_transform(P1, P2, P3, P4)
        Q1p = [Float64(SUPPORT_table[i][3]), Float64(SUPPORT_table[i][4]), 0.0]
        Q2p = [Float64(SUPPORT_table[i][5]), Float64(SUPPORT_table[i][6]), 0.0]

        stif = zeros(6)
        jadd = 2
        for j in 1:6
            sv = Float64(SUPPORT_table[i][j+4+jadd])
            sv < -1e-6 && (stif[j] = rigid)
            sv > zero_spring && (stif[j] = sv)
        end
        stif_mat = diagm(0 => stif)
        coord_flag = SUPPORT_table[i][11+jadd]
        if isa(coord_flag, String) && startswith(coord_flag, "s")
            stif_mat[1:3, 1:3] = T * stif_mat[1:3, 1:3] * T'
            stif_mat[4:6, 4:6] = T * stif_mat[4:6, 4:6] * T'
        end

        pl_elem_ids = findall(j -> el_props[j, 1] == ipl, 1:size(el_props, 1))
        for ii in pl_elem_ids
            fe_id = Int(el_props[ii, 2])
            elemnodenr = FE_dat[fe_id].nodenr
            elems_row = el_nodes[ii, 1:elemnodenr]
            nodesp = (nodes[elems_row, 1:3] .- P1') * T
            if all(nodesp[:, 1] .> Q1p[1]-tol) && all(nodesp[:, 1] .< Q2p[1]+tol) &&
               all(nodesp[:, 2] .> Q1p[2]-tol) && all(nodesp[:, 2] .< Q2p[2]+tol) &&
               all(nodesp[:, 3] .> Q1p[3]-tol) && all(nodesp[:, 3] .< Q2p[3]+tol)

                if elemnodenr == 3
                    _, P1e, P2e, P3e = ct_3node_g2e(elems_row, nodes)
                    area = (P1e[1]*P2e[2]-P2e[1]*P1e[2]-P1e[1]*P3e[2]+P3e[1]*P1e[2]+P2e[1]*P3e[2]-P3e[1]*P2e[2])/2
                else
                    _, P1e, P2e, P3e, P4e = ct_4node_g2e(elems_row, nodes)
                    area  = (P1e[1]*P2e[2]-P2e[1]*P1e[2]-P1e[1]*P3e[2]+P3e[1]*P1e[2]+P2e[1]*P3e[2]-P3e[1]*P2e[2])/2
                    area += (P1e[1]*P4e[2]-P4e[1]*P1e[2]-P1e[1]*P3e[2]+P3e[1]*P1e[2]+P4e[1]*P3e[2]-P3e[1]*P4e[2])/2
                end
                stifact = stif_mat * area / 12
                for inod in 1:elemnodenr
                    nodi = nodes[elems_row[inod], 5]
                    ndofi = (nodi - 1) * dpn
                    for jnod in 1:elemnodenr
                        nodj = nodes[elems_row[jnod], 5]
                        ndofj = (nodj - 1) * dpn
                        for idof in 1:dpn, jdof in 1:dpn
                            st = stifact[idof, jdof]
                            if abs(st) > zero_spring / 1000
                                i0 += 1
                                rows[i0] = ndofi + idof; cols[i0] = ndofj + jdof
                                vals[i0] = inod == jnod ? 2*st : st
                            end
                        end
                    end
                end
            end
        end
    end

    totaldof = nodenr * dpn
    Ke_sup = sparse(rows[1:i0], cols[1:i0], vals[1:i0], totaldof, totaldof)
    return Ke_sup
end




function cond_support(Ke, PLATE_table, CONDSUP_table, FE_dat, nodes, el_nodes, el_props, nodenr, d, dpn, uncond, tolers)
    tol = tolers.zero_dist
    zero_spring = tolers.zero_spring
    rigid = maximum(Ke) * 1000

    nonzer = nodenr * dpn^2
    rows = zeros(Int, nonzer); cols = zeros(Int, nonzer); vals = zeros(nonzer)
    i0 = 0
    nsup = length(CONDSUP_table)

    # point supports
    for i in 1:nsup
        stype = CONDSUP_table[i][1]
        if stype in ("point", "p-line", "p-surf")
            ipl = CONDSUP_table[i][2]
            P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
            P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
            _, T = plate_transform(P1, P2, P3, P4)

            if stype == "point"
                Qp = [Float64(CONDSUP_table[i][3]), Float64(CONDSUP_table[i][4]), 0.0]
                Q = T * Qp + P1
                nods = findall(j -> abs(nodes[j,1]-Q[1]) < tol && abs(nodes[j,2]-Q[2]) < tol &&
                               abs(nodes[j,3]-Q[3]) < tol && nodes[j,4] > 0, 1:size(nodes,1))
                jadd = 0
            elseif stype == "p-line"
                Q1p = [Float64(CONDSUP_table[i][3]), Float64(CONDSUP_table[i][4]), 0.0]
                Q2p = [Float64(CONDSUP_table[i][5]), Float64(CONDSUP_table[i][6]), 0.0]
                nodesp = (nodes[:, 1:3] .- P1') * T
                nods = findall(j -> abs((nodesp[j,1]-Q2p[1])*(Q1p[2]-Q2p[2]) +
                               (Q1p[1]-Q2p[1])*(Q2p[2]-nodesp[j,2])) < tol &&
                               abs(nodesp[j,3]) < tol, 1:size(nodes,1))
                jadd = 2
            else
                Q1p = [Float64(CONDSUP_table[i][3]), Float64(CONDSUP_table[i][4]), 0.0]
                Q2p = [Float64(CONDSUP_table[i][5]), Float64(CONDSUP_table[i][6]), 0.0]
                nodesp = (nodes[:, 1:3] .- P1') * T
                nods = findall(j -> nodesp[j,1] >= Q1p[1]-tol && nodesp[j,2] >= Q1p[2]-tol &&
                               nodesp[j,3] >= Q1p[3]-tol && nodesp[j,1] <= Q2p[1]+tol &&
                               nodesp[j,2] <= Q2p[2]+tol && nodesp[j,3] <= Q2p[3]+tol &&
                               nodes[j,4] > tol, 1:size(nodes,1))
                jadd = 2
            end

            stif = zeros(6)
            for j in 1:6
                sv = Float64(CONDSUP_table[i][j+4+jadd])
                sv < -1e-6 && (stif[j] = rigid)
                sv > zero_spring && (stif[j] = sv)
            end
            coord_flag = CONDSUP_table[i][11+jadd]
            dir_flag = Float64(CONDSUP_table[i][12+jadd])

            for inod in nods
                ndof = (inod - 1) * dpn
                stifact = zeros(6)
                if uncond
                    stifact = stif / 1e6
                else
                    d_node = isempty(d) ? zeros(dpn) : d[ndof+1:ndof+dpn]
                    if isa(coord_flag, String) && startswith(coord_flag, "s")
                        d_node[1:3] = T' * d_node[1:3]
                        d_node[4:dpn] = T[1:dpn-3, 1:dpn-3]' * d_node[4:dpn]
                    end
                    for jdof in 1:6
                        if jdof <= dpn && d_node[jdof] * dir_flag > 0
                            stifact[jdof] = stif[jdof]
                        end
                    end
                end
                stifact_mat = diagm(0 => stifact)
                if isa(coord_flag, String) && startswith(coord_flag, "s")
                    stifact_mat[1:3, 1:3] = T * stifact_mat[1:3, 1:3] * T'
                    stifact_mat[4:6, 4:6] = T * stifact_mat[4:6, 4:6] * T'
                end
                for idof in 1:dpn, jdof in 1:dpn
                    st = stifact_mat[idof, jdof]
                    if abs(st) > zero_spring / 1000
                        i0 += 1
                        rows[i0] = ndof + idof; cols[i0] = ndof + jdof; vals[i0] = st
                    end
                end
            end
        end
    end

    # surface supports
    for i in 1:nsup
        CONDSUP_table[i][1] != "surf" && continue
        ipl = CONDSUP_table[i][2]
        P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
        P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
        _, T = plate_transform(P1, P2, P3, P4)
        Q1p = [Float64(CONDSUP_table[i][3]), Float64(CONDSUP_table[i][4]), 0.0]
        Q2p = [Float64(CONDSUP_table[i][5]), Float64(CONDSUP_table[i][6]), 0.0]
        jadd = 2

        stif = zeros(6)
        for j in 1:6
            sv = Float64(CONDSUP_table[i][j+4+jadd])
            sv < -1e-6 && (stif[j] = rigid)
            sv > zero_spring && (stif[j] = sv)
        end
        coord_flag = CONDSUP_table[i][11+jadd]
        dir_flag = Float64(CONDSUP_table[i][12+jadd])

        pl_elem_ids = findall(j -> el_props[j, 1] == ipl, 1:size(el_props, 1))
        for ii in pl_elem_ids
            fe_id = Int(el_props[ii, 2])
            elemnodenr = FE_dat[fe_id].nodenr
            elems_row = el_nodes[ii, 1:elemnodenr]
            nodesp = (nodes[elems_row, 1:3] .- P1') * T
            if all(nodesp[:, 1] .> Q1p[1]-tol) && all(nodesp[:, 1] .< Q2p[1]+tol) &&
               all(nodesp[:, 2] .> Q1p[2]-tol) && all(nodesp[:, 2] .< Q2p[2]+tol) &&
               all(nodesp[:, 3] .> Q1p[3]-tol) && all(nodesp[:, 3] .< Q2p[3]+tol)

                if elemnodenr == 3
                    _, P1e, P2e, P3e = ct_3node_g2e(elems_row, nodes)
                    area = (P1e[1]*P2e[2]-P2e[1]*P1e[2]-P1e[1]*P3e[2]+P3e[1]*P1e[2]+P2e[1]*P3e[2]-P3e[1]*P2e[2])/2
                else
                    _, P1e, P2e, P3e, P4e = ct_4node_g2e(elems_row, nodes)
                    area  = (P1e[1]*P2e[2]-P2e[1]*P1e[2]-P1e[1]*P3e[2]+P3e[1]*P1e[2]+P2e[1]*P3e[2]-P3e[1]*P2e[2])/2
                    area += (P1e[1]*P4e[2]-P4e[1]*P1e[2]-P1e[1]*P3e[2]+P3e[1]*P1e[2]+P4e[1]*P3e[2]-P3e[1]*P4e[2])/2
                end

                stifact = zeros(6)
                if uncond
                    stifact = stif / 1000
                else
                    d_elem_mat = zeros(dpn, elemnodenr)
                    for inode in 1:elemnodenr
                        id0 = (elems_row[inode] - 1) * dpn
                        d_node = isempty(d) ? zeros(dpn) : d[id0+1:id0+dpn]
                        if isa(coord_flag, String) && startswith(coord_flag, "s")
                            d_node[1:3] = T' * d_node[1:3]
                            d_node[4:dpn] = T[1:dpn-3, 1:dpn-3]' * d_node[4:dpn]
                        end
                        d_elem_mat[:, inode] = d_node
                    end
                    d_ave = vec(mean(d_elem_mat, dims=2))
                    for jdof in 1:6
                        if jdof <= dpn && d_ave[jdof] * dir_flag > 0
                            stifact[jdof] = stif[jdof]
                        end
                    end
                end
                stifact_mat = diagm(0 => stifact)
                # coordinate flag index for surface = CONDSUP_table[i][13]
                coord_flag2 = length(CONDSUP_table[i]) >= 13 ? CONDSUP_table[i][13] : "g"
                if isa(coord_flag2, String) && startswith(coord_flag2, "s")
                    stifact_mat[1:3, 1:3] = T * stifact_mat[1:3, 1:3] * T'
                    stifact_mat[4:6, 4:6] = T * stifact_mat[4:6, 4:6] * T'
                end
                stifact_mat = stifact_mat * area / 12

                for inod in 1:elemnodenr
                    nodi = Int(nodes[elems_row[inod], 5])
                    ndofi = (nodi - 1) * dpn
                    for jnod in 1:elemnodenr
                        nodj = Int(nodes[elems_row[jnod], 5])
                        ndofj = (nodj - 1) * dpn
                        for idof in 1:dpn, jdof in 1:dpn
                            st = stifact_mat[idof, jdof]
                            if abs(st) > zero_spring / 1000
                                i0 += 1
                                rows[i0] = ndofi + idof; cols[i0] = ndofj + jdof
                                vals[i0] = inod == jnod ? 2*st : st
                            end
                        end
                    end
                end
            end
        end
    end

    totaldof = nodenr * dpn
    Ke_cond = sparse(rows[1:i0], cols[1:i0], vals[1:i0], totaldof, totaldof)
    return Ke_cond
end






"""
#########################################

STRESS

#########################################
"""


function stress_uv_3n(d_elem, P1, P2, P3, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp, ifgauss)
    xe1=P1[1]; ye1=P1[2]
    xe2=P2[1]; ye2=P2[2]
    xe3=P3[1]; ye3=P3[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)

    D = xe1*ye2 - xe2*ye1 - xe1*ye3 + xe3*ye1 + xe2*ye3 - xe3*ye2

    DB = [
        (E11*t*(ye2-ye3))/D  -(E12*t*(xe2-xe3))/D  -(E11*t*(ye1-ye3))/D  (E12*t*(xe1-xe3))/D  (E11*t*(ye1-ye2))/D  -(E12*t*(xe1-xe2))/D;
        (E21*t*(ye2-ye3))/D  -(E22*t*(xe2-xe3))/D  -(E21*t*(ye1-ye3))/D  (E22*t*(xe1-xe3))/D  (E21*t*(ye1-ye2))/D  -(E22*t*(xe1-xe2))/D;
        -(G*t*(xe2-xe3))/D    (G*t*(ye2-ye3))/D     (G*t*(xe1-xe3))/D   -(G*t*(ye1-ye3))/D   -(G*t*(xe1-xe2))/D    (G*t*(ye1-ye2))/D
    ]

    npoint = ifgauss ? nGp : 3
    str_uv = zeros(3, max(npoint, 1))
    str_uv[:, 1] = DB * d_elem / t
    for i in 2:npoint
        str_uv[:, i] = str_uv[:, 1]
    end
    return str_uv
end



function stress_uv_4n(d_elem, P1, P2, P3, P4, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp_in, ifgauss)
    xe1=P1[1]; ye1=P1[2]
    xe2=P2[1]; ye2=P2[2]
    xe3=P3[1]; ye3=P3[2]
    xe4=P4[1]; ye4=P4[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)

    if ifgauss
        nGp = round(Int, sqrt(nGp_in))
        locx, locy, wei = GL4(nGp)
    else
        nGp = 2
        locx = [-1.0  1.0;  1.0 -1.0]
        locy = [-1.0 -1.0;  1.0  1.0]
    end

    ii = 0
    str_uv = zeros(3, nGp * nGp)
    for i in 1:nGp, j in 1:nGp
        ii += 1
        x = locx[i, j]; y = locy[i, j]
        Jdet_expr = xe1*ye2-xe2*ye1-xe1*ye4+xe2*ye3-xe3*ye2+xe4*ye1+xe3*ye4-xe4*ye3
        Jdet_x    = -xe1*ye3+xe3*ye1+xe1*ye4+xe2*ye3-xe3*ye2-xe4*ye1-xe2*ye4+xe4*ye2
        Jdet_y    = -xe1*ye2+xe2*ye1+xe1*ye3-xe3*ye1-xe2*ye4+xe4*ye2+xe3*ye4-xe4*ye3
        Jdet_xy   = -xe1*ye3+xe3*ye1+xe1*ye4+xe2*ye3-xe3*ye2-xe4*ye1-xe2*ye4+xe4*ye2

        denom = Jdet_expr + x*Jdet_x + y*Jdet_y + x*y*(xe1*ye3-xe3*ye1-xe1*ye4-xe2*ye3+xe3*ye2+xe4*ye1+xe2*ye4-xe4*ye2+xe1*ye3-xe3*ye1-xe1*ye4-xe2*ye3+xe3*ye2+xe4*ye1+xe2*ye4-xe4*ye2)
        # Use the standard Jacobian formulation for 4-node element
        J11 = xe1*(y/4-1/4) - xe2*(y/4-1/4) + xe3*(y/4+1/4) - xe4*(y/4+1/4)
        J12 = ye1*(y/4-1/4) - ye2*(y/4-1/4) + ye3*(y/4+1/4) - ye4*(y/4+1/4)
        J21 = xe1*(x/4-1/4) - xe2*(x/4+1/4) + xe3*(x/4+1/4) - xe4*(x/4-1/4)
        J22 = ye1*(x/4-1/4) - ye2*(x/4+1/4) + ye3*(x/4+1/4) - ye4*(x/4-1/4)
        Jmat = [J11 J12; J21 J22]
        detJ = det(Jmat)
        Ji = inv(Jmat)
        Ji11=Ji[1,1]; Ji12=Ji[1,2]; Ji21=Ji[2,1]; Ji22=Ji[2,2]

        N1x=(y-1)/4; N2x=-(y-1)/4; N3x=(y+1)/4; N4x=-(y+1)/4
        N1y=(x-1)/4; N2y=-(x+1)/4; N3y=(x+1)/4; N4y=-(x-1)/4

        Buv = [
            N1x*Ji11+N1y*Ji12  0.0  N2x*Ji11+N2y*Ji12  0.0  N3x*Ji11+N3y*Ji12  0.0  N4x*Ji11+N4y*Ji12  0.0;
            0.0  N1x*Ji21+N1y*Ji22  0.0  N2x*Ji21+N2y*Ji22  0.0  N3x*Ji21+N3y*Ji22  0.0  N4x*Ji21+N4y*Ji22;
            N1x*Ji21+N1y*Ji22  N1x*Ji11+N1y*Ji12  N2x*Ji21+N2y*Ji22  N2x*Ji11+N2y*Ji12  N3x*Ji21+N3y*Ji22  N3x*Ji11+N3y*Ji12  N4x*Ji21+N4y*Ji22  N4x*Ji11+N4y*Ji12
        ]
        D_mat = [E11 E12 0.0; E21 E22 0.0; 0.0 0.0 G]
        str_uv[:, ii] = D_mat * Buv * d_elem / t
    end
    return str_uv
end



function main_stressout(d, PLATE_table, CONDSUP_table, MAT_table, nodes, el_nodes, el_props, FE_dat, stressout, comp, mat_sym, dpn)
    println("Calculation for visualization ")
    t0 = time()
    nodenr = length(d) ÷ dpn
    nelems = size(el_nodes, 1)
    nres = 10
    out_res = zeros(nelems, 9, nres, 3)

    ifgauss = false
    nGp = 1

    # conditional support status
    stifactive = zeros(nodenr)
    if stressout != 0 && !isempty(intersect(comp, [10]))
        nsup = length(CONDSUP_table)
        tol = 1e-4
        for i in 1:nsup
            CONDSUP_table[i][1] != "surf" && continue
            ipl = CONDSUP_table[i][2]
            Q1p = [Float64(CONDSUP_table[i][3]), Float64(CONDSUP_table[i][4]), 0.0]
            Q2p = [Float64(CONDSUP_table[i][5]), Float64(CONDSUP_table[i][6]), 0.0]
            P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
            P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
            _, T = plate_transform(P1, P2, P3, P4)
            nodesp = (nodes[:, 1:3] .- P1') * T
            nods = findall(j -> nodesp[j,1] >= Q1p[1]-tol && nodesp[j,2] >= Q1p[2]-tol &&
                           nodesp[j,3] >= Q1p[3]-tol && nodesp[j,1] <= Q2p[1]+tol &&
                           nodesp[j,2] <= Q2p[2]+tol && nodesp[j,3] <= Q2p[3]+tol &&
                           nodes[j,4] > tol, 1:size(nodes,1))
            dir_flag = Float64(CONDSUP_table[i][14])
            for inod in nods
                nod = inod
                ndof = (nod - 1) * dpn
                if d[ndof+3] * dir_flag > 0
                    stifactive[nod] = dir_flag
                end
            end
        end
        for ii in 1:nelems
            fe_id = Int(el_props[ii, 2]); fe_id == 0 && continue
            elemnodenr = FE_dat[fe_id].nodenr
            for j in 1:elemnodenr
                n0 = el_nodes[ii, j]
                out_res[ii, j, 10, 2] = stifactive[n0]
            end
        end
    end

    # stresses
    if stressout != 0 && !isempty(intersect(comp, 7:9))
        for ii in 1:nelems
            fe_id = Int(el_props[ii, 2]); fe_id == 0 && continue
            elemnodenr = FE_dat[fe_id].nodenr
            ipl = Int(el_props[ii, 1])
            P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
            P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
            _, Tgs = plate_transform(P1, P2, P3, P4)
            tel = el_props[ii, 3]; imat = Int(el_props[ii, 4])
            mv = MAT_table[imat]
            Ex=mv[1]; Ey=mv[2]; nuxy=mv[3]; nuyx=mv[4]; G=mv[5]
            elems_row = el_nodes[ii, 1:elemnodenr]

            if elemnodenr == 3
                Tge, P1e, P2e, P3e = ct_3node_g2e(elems_row, nodes)
            else
                Tge, P1e, P2e, P3e, P4e = ct_4node_g2e(elems_row, nodes)
            end

            d_elem = Float64[]; induv = Int[]
            for inode in 1:elemnodenr
                i0_n = (elems_row[inode] - 1) * dpn
                d_node = copy(d[i0_n+1:i0_n+dpn])
                d_node[1:3] = Tge' * d_node[1:3]
                d_node[4:dpn] = Tge[1:dpn-3, 1:dpn-3]' * d_node[4:dpn]
                append!(d_elem, d_node)
                append!(induv, [(inode-1)*dpn+1, (inode-1)*dpn+2])
            end

            if elemnodenr == 3
                str_uv = stress_uv_3n(d_elem[induv], P1e, P2e, P3e, tel, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp, ifgauss)
            else
                str_uv = stress_uv_4n(d_elem[induv], P1e, P2e, P3e, P4e, tel, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp, ifgauss)
            end

            Tse = Tgs' * Tge
            c = Tse[1,1]; s = Tse[2,1]
            Tstress = [c^2 s^2 2*s*c; s^2 c^2 -2*s*c; -s*c s*c c^2-s^2]
            str_uv = Tstress * str_uv

            for j in 1:elemnodenr
                out_res[ii, j, 7:9, 2] = str_uv[1:3, j]
            end
        end
    end

    # displacements
    if stressout != 0 && !isempty(intersect(comp, 1:6))
        for ii in 1:nelems
            fe_id = Int(el_props[ii, 2]); fe_id == 0 && continue
            elemnodenr = FE_dat[fe_id].nodenr
            ipl = Int(el_props[ii, 1])
            P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
            P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
            _, Tgs = plate_transform(P1, P2, P3, P4)
            for j in 1:elemnodenr
                n0 = el_nodes[ii, j]
                i0_n = (n0 - 1) * dpn
                dXYZ = d[i0_n+1:i0_n+3]
                out_res[ii, j, 1:3, 2] = dXYZ
                out_res[ii, j, 4:6, 2] = Tgs' * dXYZ
            end
        end
    end

    println("  done in $(round(time()-t0, digits=2)) s")
    return out_res
end




"""
#########################################

ANALYSIS

#########################################
"""


function FEanal_incr(iincr, dincr, d_ini, Ke, lv, nodes, el_nodes, el_props,
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
     #   lv = load_vec_UDL(nodes, el_nodes, el_props, nodenr, PLATE_table, LOAD_table, FE_dat, dpn)
     #   lv .+= load_vec_edge(nodes, el_nodes, el_props, pl_el_edges, pl_edge_nodes, nodenr, PLATE_table, LOAD_table, FE_dat, dpn, tolers)
     #   lv .+= load_vec_point(nodes, nodenr, PLATE_table, LOAD_table, dpn, tolers)
        lv = load_vec_edge(nodes, el_nodes, el_props, pl_el_edges, pl_edge_nodes, nodenr, PLATE_table, LOAD_table, FE_dat, dpn, tolers)
        println("  done in $(round(time()-t0, digits=2)) s")
    end

    println("Load increment $iincr")
    println("Cond supports")
    t0 = time()
    dlv = lv * dincr[iincr]

    if iincr == 1
        uncond = true
        Ke_cond = cond_support(Ke, PLATE_table, CONDSUP_table, FE_dat, nodes, el_nodes, el_props, nodenr, d_ini, dpn, uncond, tolers)
        d_ini = Symmetric(Ke + Ke_cond) \ Array(dlv)
    end

    uncond = false
    Ke_cond = cond_support(Ke, PLATE_table, CONDSUP_table, FE_dat, nodes, el_nodes, el_props, nodenr, d_ini, dpn, uncond, tolers)
    println("  done in $(round(time()-t0, digits=2)) s")
    t0 = time()
    println("Equation solve")
    Kesolve=Ke + Ke_cond
   # KeT=Kesolve'
   # Kesolve=KeT'
 #   dd = Array(Ke + Ke_cond) \ dlv
    dd = Symmetric(Kesolve) \ Array(dlv)
    println("  done in $(round(time()-t0, digits=2)) s")
    if iincr == 1; d_ini = zeros(length(d_ini)); end
    d_stat = d_ini + dd

    return d_stat, Ke, lv
end



"""
#########################################

VISUALIZATION, OUTPUT

#########################################
"""


function plot_res_tri(el_nodes, el_props, out_grid, FE_dat, out_res, scale, layer, comp, gridon;
                      title_str = "Component $comp",
                      colormap  = :jet,
                      azimuth   = deg2rad(60),
                      elevation = deg2rad(20))

    nel    = size(el_nodes, 1)
    nnodes = size(out_grid, 1)

    # ── single pass: per-node displacement + result value ────────────────────
    node_disp = zeros(nnodes, 3)
    node_val  = zeros(nnodes)
    node_cnt  = zeros(Int, nnodes)
    node_seen = falses(nnodes)

    tri_list = NTuple{3,Int}[]     # triangulated face connectivity

    for iel in 1:nel
        fe_id = Int(el_props[iel, 2])
        fe_id == 0 && continue
        elemnodenr = FE_dat[fe_id].nodenr

        for j in 1:elemnodenr
            n = el_nodes[iel, j]
            n == 0 && continue
            if !node_seen[n]
                node_disp[n, :] = out_res[iel, j, 1:3, 2]
                node_seen[n] = true
            end
            node_val[n] += out_res[iel, j, comp, layer]
            node_cnt[n] += 1
        end

        el_props[iel, 5] > 0 || continue        # skip hole elements

        if elemnodenr == 3
            push!(tri_list, (el_nodes[iel,1], el_nodes[iel,2], el_nodes[iel,3]))
        elseif elemnodenr == 4
            push!(tri_list, (el_nodes[iel,1], el_nodes[iel,2], el_nodes[iel,3]))
            push!(tri_list, (el_nodes[iel,1], el_nodes[iel,3], el_nodes[iel,4]))
        end
    end

    # ── per-node averaged result ──────────────────────────────────────────────
    out_comp = [node_cnt[i] > 0 ? node_val[i] / node_cnt[i] : 0.0 for i in 1:nnodes]

    # ── displaced positions ───────────────────────────────────────────────────
    def_grid = out_grid[:, 1:3] .+ scale .* node_disp

    # ── GeometryBasics mesh (required by both mesh! and wireframe!) ───────────
    verts = [Point3f(def_grid[i, 1], def_grid[i, 2], def_grid[i, 3]) for i in 1:nnodes]
    faces = [TriangleFace(t[1], t[2], t[3]) for t in tri_list]
    geom  = Mesh(verts, faces)

#     # ── figure ────────────────────────────────────────────────────────────────
#     fig = Figure(size = (800, 600))
#     ax  = Axis3(fig[1, 1];
#                 aspect    = :data,
#                 xlabel    = "X (in)",
#                 ylabel    = "Y (in)",
#                 zlabel    = "Z (in)",
#                 title     = title_str,
#                 azimuth   = azimuth,
#                 elevation = elevation)

#     # filled colour-mapped mesh with smooth per-vertex interpolation
#     msh = mesh!(ax, geom;
#                 color    = out_comp,
#                 colormap = colormap,
#                 shading  = NoShading)

#     # optional grey wireframe edge overlay
#     if gridon
#         wireframe!(ax, geom;
#                    color     = (:gray, 0.35),
#                    linewidth = 0.5)
#     end

#     Colorbar(fig[1, 2], msh; width = 20, label = title_str)

#     return fig
# end













end  # module 

using BasePlateFEA
using SparseArrays



# ─────────────────────────────────────────────────────────────────────────────
# SETTINGS
# ─────────────────────────────────────────────────────────────────────────────
dpn = 6
calcmodenr = 20
stressout = 0

# tolerances (kip-inch units after conversion)
tolers = BasePlateFEA.Tolerances(
    1e-3 / 25.4,          # zero_thickn
    1e-3 / 4448.222 * 25.4,  # zero_spring
    1e-2 / 25.4,          # zero_dist
    0.2  / 25.4            # zeroFElength
)




# finite element properties
FE_dat = BasePlateFEA.FE_props()

# ─────────────────────────────────────────────────────────────────────────────
# GEOMETRY: C-member with base-plate
# ─────────────────────────────────────────────────────────────────────────────
W = 5.0; L = 5.5; t = 0.375; Da = 0.75; Dc = 0.25; Sa = 3.5
R = 1/4  # bolt-hole radius

h = 3.0; b = h; c = 0.5
tb = 0.0747; th = tb; tc = tb
H = 3.0

A_cs, yCG, zCG, Iyr, Iy, Izr, Iz = BasePlateFEA.cs_data(h, b, c, tb, th, tc)

# ─────────────────────────────────────────────────────────────────────────────
# MESH control
# ─────────────────────────────────────────────────────────────────────────────
fe_size = [0.15, 0.15 ]
nplate = 6
mesh_typ = fill(0, nplate, 2)
for i in 1:nplate
    mesh_typ[i, :] = [12, 22]
end

rat_q2t = [0.0, 1.0, 1.0, 1.0, 1.0, 1.0]

# ─────────────────────────────────────────────────────────────────────────────
# MATERIAL
# ─────────────────────────────────────────────────────────────────────────────
E = 29500.0; nu = 0.3; G = E / 2 / (1 + nu); ro = 1.0
MAT_table = [
    [E, E, nu, nu, G, ro],
    [E, E, nu, nu, G, ro],
]
mat_sym = 1

# ─────────────────────────────────────────────────────────────────────────────
# PLATE TABLE  [x1 y1 z1  x2 y2 z2  x3 y3 z3  x4 y4 z4  t  imat]
# ─────────────────────────────────────────────────────────────────────────────
PLATE_table = [
    [0.0, -W/2, 0.0,  L,   -W/2, 0.0,  L,    W/2,    0.0,  0.0,  W/2,    0.0,  t,   1.0],
    [Dc+b, h/2-c, 0.0, Dc+b, h/2, 0.0, Dc+b, h/2,    H,    Dc+b, h/2-c,  H,    tc,  1.0],
    [Dc+b, h/2, 0.0,  Dc,   h/2, 0.0,  Dc,   h/2,    H,    Dc+b, h/2,    H,    tb,  1.0],
    [Dc,   h/2, 0.0,  Dc,  -h/2, 0.0,  Dc,  -h/2,    H,    Dc,   h/2,    H,    th,  1.0],
    [Dc,  -h/2, 0.0,  Dc+b,-h/2, 0.0,  Dc+b,-h/2,    H,    Dc,  -h/2,   H,    tb,  1.0],
    [Dc+b,-h/2, 0.0,  Dc+b,-h/2+c,0.0, Dc+b,-h/2+c,  H,    Dc+b,-h/2,   H,    tc,  1.0],
]

# shift y coordinates so plate 1 sits at y=0..W
for ipl in 1:length(PLATE_table)
    PLATE_table[ipl][2:3:11] .+= W/2
end

# ─────────────────────────────────────────────────────────────────────────────
# LOADS
# ─────────────────────────────────────────────────────────────────────────────
N0 = -1.0; N = 0.5 * N0; p_cs = N / A_cs

LOAD_table = Vector{Vector{Any}}(undef, 18)
# normal force (axial)
i0 = 0
LOAD_table[i0+1] = Any["edge", 2, 0.0, H, c, H, 0.0, 0.0, p_cs*tc, 0.0, 0.0, p_cs*tc, "gl"]
LOAD_table[i0+2] = Any["edge", 3, 0.0, H, b, H, 0.0, 0.0, p_cs*tc, 0.0, 0.0, p_cs*tb, "gl"]
LOAD_table[i0+3] = Any["edge", 4, 0.0, H, h, H, 0.0, 0.0, p_cs*tc, 0.0, 0.0, p_cs*th, "gl"]
LOAD_table[i0+4] = Any["edge", 5, 0.0, H, b, H, 0.0, 0.0, p_cs*tc, 0.0, 0.0, p_cs*tb, "gl"]
LOAD_table[i0+5] = Any["edge", 6, 0.0, H, c, H, 0.0, 0.0, p_cs*tc, 0.0, 0.0, p_cs*tc, "gl"]

# horizontal force along symmetry axis
V2 = 0.0; p_v = V2 / b
i0 = 5
LOAD_table[i0+1] = Any["edge", 3, 0.0, H, b, H, p_v, 0.0, 0.0, p_v, 0.0, 0.0, "gl"]
LOAD_table[i0+2] = Any["edge", 5, 0.0, H, b, H, p_v, 0.0, 0.0, p_v, 0.0, 0.0, "gl"]

# horizontal force perpendicular to symmetry axis
V1 = 0.0; p_v1 = V1 / h
i0 = 7
LOAD_table[i0+1] = Any["edge", 4, 0.0, H, h, H, 0.0, p_v1, 0.0, 0.0, p_v1, 0.0, "gl"]

# major axis bending
M0 = 0.0; M1 = -1.0 * M0
p1b = M1 / Iy * (yCG - b); p2b = M1 / Iy * yCG
i0 = 8
LOAD_table[i0+1] = Any["edge", 2, 0.0, H, c, H, 0.0, 0.0, p1b*tc, 0.0, 0.0, p1b*tc, "gl"]
LOAD_table[i0+2] = Any["edge", 3, 0.0, H, b, H, 0.0, 0.0, p1b*tb, 0.0, 0.0, p2b*tb, "gl"]
LOAD_table[i0+3] = Any["edge", 4, 0.0, H, h, H, 0.0, 0.0, p2b*th, 0.0, 0.0, p2b*th, "gl"]
LOAD_table[i0+4] = Any["edge", 5, 0.0, H, b, H, 0.0, 0.0, p2b*tb, 0.0, 0.0, p1b*tb, "gl"]
LOAD_table[i0+5] = Any["edge", 6, 0.0, H, c, H, 0.0, 0.0, p1b*tc, 0.0, 0.0, p1b*tc, "gl"]

# minor axis bending
M0 = 1.0; M2 = 1.0 * M0
p1m = M2 / Iz * h/2; p2m = M2 / Iz * (h/2 - c)
i0 = 13
LOAD_table[i0+1] = Any["edge", 2, 0.0, H, c, H, 0.0, 0.0, p2m*tc, 0.0, 0.0, p1m*tc, "gl"]
LOAD_table[i0+2] = Any["edge", 3, 0.0, H, b, H, 0.0, 0.0, p1m*tb, 0.0, 0.0, p1m*tb, "gl"]
LOAD_table[i0+3] = Any["edge", 4, 0.0, H, h, H, 0.0, 0.0, p1m*th, 0.0, 0.0, -p1m*th, "gl"]
LOAD_table[i0+4] = Any["edge", 5, 0.0, H, b, H, 0.0, 0.0, -p1m*tb, 0.0, 0.0, -p1m*tb, "gl"]
LOAD_table[i0+5] = Any["edge", 6, 0.0, H, c, H, 0.0, 0.0, -p1m*tc, 0.0, 0.0, -p2m*tc, "gl"]

# ─────────────────────────────────────────────────────────────────────────────
# BOLT HOLES
# ─────────────────────────────────────────────────────────────────────────────
HOLE_table = Vector{Union{Nothing, BasePlateFEA.HoleData}}(nothing, nplate)

x1 = L - Da; y1 = W/2 - Sa/2
x2 = L - Da; y2 = W/2 + Sa/2
sec = 4*4  # polygon sides (multiple of 4)
dalp = 2*pi / sec
alp0 = 0.0

holes1 = zeros(sec+1, 2); holes2 = zeros(sec+1, 2)
for i in 1:sec+1
    alp = alp0 + (i-1)*dalp
    xp = R*cos(alp); yp = R*sin(alp)
    holes1[i, :] = [x1-xp, y1-yp]
    holes2[i, :] = [x2-xp, y2-yp]
end
HOLE_table[1] = BasePlateFEA.HoleData(
    [[x1, y1], [x2, y2]],
    [R, R],
    [sec, sec],
    [holes1, holes2]
)

# ─────────────────────────────────────────────────────────────────────────────
# SUPPORTS
# ─────────────────────────────────────────────────────────────────────────────
rx = 100*dalp*R; ry = 100*dalp*R; rz_bolt = 100000*dalp*R

SUPPORT_table  = Vector{Vector{Any}}()
CONDSUP_table  = Vector{Vector{Any}}()

i0 = 0
for i in 1:sec
    alp = alp0 + (i-1)*dalp
    xp = R*cos(alp); yp = R*sin(alp)
    push!(CONDSUP_table, Any["point", 1, x1-xp, y1-yp, 0.0, 0.0, rz_bolt, 0.0, 0.0, 0.0, "g", 1.0])
    push!(CONDSUP_table, Any["point", 1, x2-xp, y2-yp, 0.0, 0.0, rz_bolt, 0.0, 0.0, 0.0, "g", 1.0])
    push!(SUPPORT_table,  Any["point", 1, x1-xp, y1-yp, rx, ry, 0.0, 0.0, 0.0, 0.0, "g"])
    push!(SUPPORT_table,  Any["point", 1, x2-xp, y2-yp, ry, ry, 0.0, 0.0, 0.0, 0.0, "g"])
end

# under-baseplate conditional support (compression only)
push!(CONDSUP_table, Any["surf", 1, 0.0, 0.0, L, W, 0.0, 0.0, 10000.0, 0.0, 0.0, 0.0, "g", -1.0])

# ─────────────────────────────────────────────────────────────────────────────
# ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────
nodes, el_nodes, el_props, nodenr, pl_el_edges, pl_edge_nodes =
    BasePlateFEA.main_modelgen(PLATE_table, LOAD_table, SUPPORT_table, CONDSUP_table, HOLE_table,
                  FE_dat, fe_size, mesh_typ, rat_q2t, tolers)

d_act = Float64[]; lv = Float64[]; Ke = spzeros(0, 0); Ke_cond = spzeros(0, 0);




# incremental load steps
# lf = vcat(collect(0.000:0.001:0.005), collect(0.01:0.01:0.05), collect(0.1:0.1:1.0))
# lf = collect(0.000:0.001:0.004)
lf = [0, 0.001, 0.002, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1]
dincr = diff(lf)
nincr = length(dincr)
dmami = zeros(nincr+1, 2)
d_act=zeros(nodenr*dpn)

for iincr in 1:nincr
#    global d_ini, Ke, lv
#    d_stat, Ke, lv = BasePlateFEA.FEanal_incr(iincr, dincr, d_ini, Ke, lv,
#                                   nodes, el_nodes, el_props, pl_el_edges, pl_edge_nodes,
#                                   nodenr, PLATE_table, MAT_table, SUPPORT_table,
#                                   CONDSUP_table, LOAD_table, FE_dat, mat_sym, dpn, tolers)
#    d_ini = d_stat   # keep at nodenr*dpn for next FEanal_incr call

   global d_act, Ke, Ke_cond, lv
   d_act, Ke, Ke_cond, lv = BasePlateFEA.FEanal_incr(iincr, dincr, d_act, Ke, Ke_cond, lv,
                                   nodes, el_nodes, el_props, pl_el_edges, pl_edge_nodes,
                                   nodenr, PLATE_table, MAT_table, SUPPORT_table,
                                   CONDSUP_table, LOAD_table, FE_dat, mat_sym, dpn, tolers)

    # re-index to all-node DOF space for post-processing
    ind = Int.(nodes[:, 5]) .- 1
    ndof_full = length(ind) * dpn
    ind2 = zeros(Int, ndof_full)
    for i in 1:dpn
        ind2[i:dpn:ndof_full] = ind .* dpn .+ i
    end
    d_stat = d_act[ind2]   # now total_nodes*dpn (for dmami, visualization)

    ind_bp = findall(j -> abs(nodes[j, 3]) < 0.001, 1:size(nodes, 1))
    dmami[iincr+1, 1] = minimum(d_stat[(ind_bp .- 1) .* dpn .+ 3])
    dmami[iincr+1, 2] = maximum(d_stat[(ind_bp .- 1) .* dpn .+ 3])

    # visualization at last increment
    if stressout > 0.01 && iincr == nincr
        d_out = d_stat
        comp_vis = collect(1:10)
        out_res = BasePlateFEA.main_stressout(d_out, PLATE_table, CONDSUP_table, MAT_table,
                                  nodes, el_nodes, el_props, FE_dat, stressout,
                                  comp_vis, mat_sym, dpn)

        scale_vis = max(W, L) / max(maximum(abs.(d_out)), 1e-30) / 10
        layer = 2

        # println("Plotting displacement (comp=3)...")
        # fig1 = BasePlateFEA.plot_res_tri(el_nodes, el_props, nodes, FE_dat, out_res, scale_vis, layer, 3, true;
        #                     title_str="Vertical Displacement W (in)")
        # display(fig1)
        # save(joinpath(@__DIR__, "..", "displacement_w.png"), fig1)

      #  println("Plotting stress σ_y (comp=8)...")
      #  fig2 = plot_res_tri(el_nodes, el_props, nodes, FE_dat, out_res, scale_vis, layer, 8, true;
      #                      title_str="In-Plane Stress σ_y (ksi)")
      #  save(joinpath(@__DIR__, "..", "stress_sy.png"), fig2)

      #  println("Plotting contact (comp=10)...")
      #  fig3 = plot_res_tri(el_nodes, el_props, nodes, FE_dat, out_res, 0.0, layer, 10, true;
      #                      title_str="Contact Status")
      #  save(joinpath(@__DIR__, "..", "contact.png"), fig3)
    end
end

# # ─────────────────────────────────────────────────────────────────────────────
# # FORCE-DISPLACEMENT CURVE
# # ─────────────────────────────────────────────────────────────────────────────
zerodist = tolers.zero_dist
if abs(M2) > 0.00001
    nod_b1 = findall(j -> abs(nodes[j,3]) < zerodist &&
                     nodes[j,1] > Dc+yCG-fe_size[1] && nodes[j,1] < Dc+yCG+fe_size[1] &&
                     abs(nodes[j,2] - (W/2 - h/2)) < zerodist && nodes[j,4] > 0,
                     1:size(nodes,1))
    nod_b2 = findall(j -> abs(nodes[j,3]) < zerodist &&
                     nodes[j,1] > Dc+yCG-fe_size[1] && nodes[j,1] < Dc+yCG+fe_size[1] &&
                     abs(nodes[j,2] - (W/2 + h/2)) < zerodist && nodes[j,4] > 0,
                     1:size(nodes,1))
    if !isempty(nod_b1) && !isempty(nod_b2)
        act_b1 = Int.(nodes[nod_b1, 5])
        act_b2 = Int.(nodes[nod_b2, 5])
        dof_b1 = (d_act[(act_b1 .- 1) .* dpn .+ 3])
        dof_b2 = (d_act[(act_b2 .- 1) .* dpn .+ 3])
        ave_d_b1=sum(dof_b1)/size(dof_b1,1)
        ave_d_b2=sum(dof_b2)/size(dof_b2,1)
        rot  = (ave_d_b2 - ave_d_b1) / h
        stif2 = abs(M2) / abs(rot)
        println("Rotational stiffness = $stif2 kip-in/rad")
    end
end

# nod_bp = findall(j -> abs(nodes[j, 3]) < zerodist, 1:size(nodes, 1))
# act_bp = Int.(nodes[nod_bp, 5])
# d_bp_max = maximum(d_ini[(act_bp .- 1) .* dpn .+ 3])



# # save force-displacement data as CSV
# using DelimitedFiles
# csv_path = joinpath(@__DIR__, "..", "force_disp_data.csv")
# header = ["load_factor" "disp_min_in" "disp_max_in" "moment_lbin"]
# data   = hcat(lf, dmami[:, 1], dmami[:, 2], M2 * 1000 .* lf)
# writedlm(csv_path, [header; data], ',')
# println("Saved force_disp_data.csv")

# println("End")
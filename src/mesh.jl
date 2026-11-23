# --------------------------------------------------------------------------------------------
# Model generation (main_modelgen_dm.m, mesh4_rect.m, collect_*.m): quadrilateral elements on the
# upright walls, DistMesh triangles on the base plate with the wall footprint, the hole polygons
# and the hole centres as fixed nodes, soft elements inside the holes, and the rigid-link beams
# (top cross-section to a master node, hole centre to hole perimeter).
# --------------------------------------------------------------------------------------------

"""
    Mesh

- `x`: `nnodes × 3` node coordinates
- `conn`: element connectivity (3 or 4 node numbers, counter-clockwise)
- `part`: plate number of each element (1 = base plate; 2…6 = lip, flange, web, flange, lip)
- `thickness`, `hole`: element thickness (already reduced inside the holes) and hole marker
- `beams`: `(node1, node2, kind)` rigid links, `kind` 1 = top cross-section, 2 = hole
- `bp_nodes`: base plate nodes (each has a concrete dof and a contact link, in this order)
- `bp_area`: tributary area of the base plate nodes
- `master`: master node of the top cross-section (loads are applied here)
- `bolt_nodes`: hole centre nodes (anchor bolt springs); `hole_nodes`: perimeter nodes of each hole
- `rot_nodes`: the two flange nodes between which the rotation is measured
- `corner_nodes`: base plate corner nodes `(0,0)`, `(L,0)`, `(L,W)`, `(0,W)`
- `frames`: 3×3 plate frames `[j1 j2 j3]` (columns) used to report stresses, one per part
"""
struct Mesh
    x::Matrix{Float64}
    conn::Vector{Vector{Int}}
    part::Vector{Int}
    thickness::Vector{Float64}
    hole::BitVector
    beams::Vector{NTuple{3,Int}}
    bp_nodes::Vector{Int}
    bp_area::Vector{Float64}
    master::Int
    bolt_nodes::Vector{Int}
    hole_nodes::Vector{Vector{Int}}
    rot_nodes::NTuple{2,Int}
    corner_nodes::NTuple{4,Int}
    frames::Vector{Matrix{Float64}}
end

nnodes(m::Mesh) = size(m.x, 1)
nelements(m::Mesh) = length(m.conn)

# seeds along a segment of length `len` with interior anchor positions (mesh4_rect.m)
function _seeds(len, anchors, size, tol)
    pts = [0.0, len]
    for a in anchors
        (a > tol && a < len - tol && minimum(abs.(pts .- a)) > tol) && push!(pts, a)
    end
    sort!(pts)
    s = [pts[1]]
    for i in 2:length(pts)
        n = max(round(Int, (pts[i] - pts[i-1]) / size), 1)
        d = (pts[i] - pts[i-1]) / n
        for k in 1:n
            push!(s, pts[i-1] + k * d)
        end
    end
    return s
end

"plan-view start and end points of the five upright walls (plates 2…6 of `bp_Barcelona.m`) and their extra seed positions"
function upright_walls(g::BasePlateGeometry, ym)
    h, b, c = centerline(g)
    yc = g.W / 2
    Dc = g.Dc
    return [
        (p1 = (Dc + b, yc + h / 2 - c), p2 = (Dc + b, yc + h / 2), anchors = Float64[]),
        (p1 = (Dc + b, yc + h / 2), p2 = (Dc, yc + h / 2), anchors = [b - ym]),
        (p1 = (Dc, yc + h / 2), p2 = (Dc, yc - h / 2), anchors = [c, h - c]),
        (p1 = (Dc, yc - h / 2), p2 = (Dc + b, yc - h / 2), anchors = [ym]),
        (p1 = (Dc + b, yc - h / 2), p2 = (Dc + b, yc - h / 2 + c), anchors = Float64[]),
    ]
end

"centres of the two bolt holes"
hole_centers(g::BasePlateGeometry) = [(g.L - g.Da, g.W / 2 - g.Sa / 2), (g.L - g.Da, g.W / 2 + g.Sa / 2)]

"polygon vertices of the hole with centre `c` (as in `bp_Barcelona.m`, starting on the −x side)"
hole_polygon(g::BasePlateGeometry, c) = [(c[1] - g.R * cos(2π * (i - 1) / g.hole_sectors), c[2] - g.R * sin(2π * (i - 1) / g.hole_sectors))
                                         for i in 1:g.hole_sectors]

"""
    generate_mesh(inp::BasePlateInput) -> Mesh
"""
function generate_mesh(inp::BasePlateInput)
    g = inp.geometry
    a = inp.fe_size
    tol = 1e-2 / 25.4                       # tolers.zero_dist
    H, _, Hm = rounded_heights(inp)
    h, b, c = centerline(g)
    ym = inp.ym_ratio * b

    xs = NTuple{3,Float64}[]                # nodes, merged by coordinates (within tol)
    cells = Dict{NTuple{3,Int},Vector{Int}}()
    cellof(q) = (floor(Int, q[1] / tol), floor(Int, q[2] / tol), floor(Int, q[3] / tol))
    function find(q)
        c = cellof(q)
        for dx in -1:1, dy in -1:1, dz in -1:1
            for i in get(cells, (c[1] + dx, c[2] + dy, c[3] + dz), ())
                (abs(xs[i][1] - q[1]) < tol && abs(xs[i][2] - q[2]) < tol && abs(xs[i][3] - q[3]) < tol) && return i
            end
        end
        return 0
    end
    function node!(q)
        i = find(q)
        i > 0 && return i
        push!(xs, q)
        push!(get!(cells, cellof(q), Int[]), length(xs))
        return length(xs)
    end

    conn = Vector{Int}[]; part = Int[]; thickness = Float64[]; frames = Matrix{Float64}[]
    push!(frames, Matrix(1.0I, 3, 3))

    # ---- upright walls: structured quads -------------------------------------------------
    if inp.H_fine < H - tol                 # fine up to H_fine, element height fe_size_upper above
        Hf = max(round(inp.H_fine / a) * a, Hm)
        zs = _seeds(Hf, Float64[], a, tol)
        append!(zs, Hf .+ _seeds(H - Hf, Float64[], inp.fe_size_upper, tol)[2:end])
    else
        zs = _seeds(H, Float64[], a, tol)
    end
    footprint = NTuple{2,Float64}[]
    walls = upright_walls(g, ym)
    wall_xy = Vector{Vector{NTuple{2,Float64}}}()
    for w in walls
        len = hypot(w.p2[1] - w.p1[1], w.p2[2] - w.p1[2])
        e = ((w.p2[1] - w.p1[1]) / len, (w.p2[2] - w.p1[2]) / len)
        ss = _seeds(len, w.anchors, a, tol)
        push!(wall_xy, [(w.p1[1] + s * e[1], w.p1[2] + s * e[2]) for s in ss])
        append!(footprint, wall_xy[end])
        j1 = [e[1], e[2], 0.0]; j2 = [0.0, 0.0, 1.0]
        push!(frames, hcat(j1, j2, cross(j1, j2)))
    end

    # ---- base plate: DistMesh triangles ---------------------------------------------------
    centers = hole_centers(g)
    polys = [hole_polygon(g, cc) for cc in centers]
    pfix = vcat([(0.0, 0.0), (g.L, 0.0), (g.L, g.W), (0.0, g.W)], footprint, polys..., centers)
    fd = q -> -min(q[2], g.W - q[2], q[1], g.L - q[1])
    p, t = distmesh2d(fd, a, ((0.0, 0.0), (g.L, g.W)), pfix)
    pid = [node!((q[1], q[2], 0.0)) for q in p]
    Rin = g.R * cos(π / g.hole_sectors)
    hole = falses(0)
    for T in t
        push!(conn, [pid[T[1]], pid[T[2]], pid[T[3]]]); push!(part, 1)
        cx = (p[T[1]][1] + p[T[2]][1] + p[T[3]][1]) / 3; cy = (p[T[1]][2] + p[T[2]][2] + p[T[3]][2]) / 3
        inhole = any(hypot(cx - cc[1], cy - cc[2]) < Rin + tol for cc in centers)
        push!(hole, inhole)
        push!(thickness, inhole ? g.t * inp.hole_thickness_factor : g.t)
    end
    nbp = length(xs)                         # the base plate nodes are numbered first

    for (iw, xy) in enumerate(wall_xy)
        ids = [node!((q[1], q[2], z)) for q in xy, z in zs]
        for j in 1:length(zs)-1, i in 1:length(xy)-1
            push!(conn, [ids[i, j], ids[i+1, j], ids[i+1, j+1], ids[i, j+1]])
            push!(part, iw + 1); push!(thickness, g.tu); push!(hole, false)
        end
    end
    length(xs) > nbp && all(q[3] > tol for q in xs[nbp+1:end]) ||
        error("the upright footprint nodes were not all found in the base plate mesh")

    function findnode(q, what)
        i = find(q)
        i == 0 && error("mesh generation: no node at $what $q")
        return i
    end

    # ---- rigid links --------------------------------------------------------------------
    beams = NTuple{3,Int}[]
    top = [i for i in 1:length(xs) if abs(xs[i][3] - H) < tol]
    nprev = length(xs)
    master = node!((g.Dc + b / 2, g.W / 2, H))
    master == nprev + 1 || error("the master node coincides with a mesh node")
    for i in top
        push!(beams, (master, i, 1))
    end
    bolt_nodes = [findnode((cc[1], cc[2], 0.0), "hole centre") for cc in centers]
    hole_nodes = [[findnode((q[1], q[2], 0.0), "hole perimeter point") for q in poly] for poly in polys]
    for (ih, m) in enumerate(bolt_nodes), i in hole_nodes[ih]
        push!(beams, (m, i, 2))
    end

    # ---- control nodes ------------------------------------------------------------------
    yc = g.W / 2
    rot_nodes = (findnode((g.Dc + ym, yc - h / 2, Hm), "rotation measurement point"),
                 findnode((g.Dc + ym, yc + h / 2, Hm), "rotation measurement point"))
    corners = (findnode((0.0, 0.0, 0.0), "corner"), findnode((g.L, 0.0, 0.0), "corner"),
               findnode((g.L, g.W, 0.0), "corner"), findnode((0.0, g.W, 0.0), "corner"))

    x = [xs[i][k] for i in 1:length(xs), k in 1:3]
    bp_nodes = collect(1:nbp)
    return Mesh(x, conn, part, thickness, hole, beams, bp_nodes, nodal_areas(x, conn, part, bp_nodes), master,
                bolt_nodes, hole_nodes, rot_nodes, corners, frames)
end

"tributary areas of the base plate nodes: one third of the adjacent base plate triangles (calc_bp_area.m)"
function nodal_areas(x, conn, part, bp_nodes)
    area = zeros(size(x, 1))
    for (e, c) in enumerate(conn)
        part[e] == 1 || continue
        A = abs((x[c[2], 1] - x[c[1], 1]) * (x[c[3], 2] - x[c[1], 2]) - (x[c[2], 2] - x[c[1], 2]) * (x[c[3], 1] - x[c[1], 1])) / 2
        for i in c
            area[i] += A / 3
        end
    end
    return area[bp_nodes]
end

"""
    mesh_from_matlab(dir, tag, inp) -> Mesh

Build the mesh from the arrays written by the Octave reference runner (`<tag>_nodes.txt`,
`_el_nodes.txt`, `_el_props.txt`, `_beam_nodes.txt`, `_beam_props.txt`, `_nods_bp.txt`,
`_nods_area.txt`, `_upright_points.txt`, `_bp_corner.txt`, `_nods_spr_b.txt`), so that the Julia
and MATLAB/Octave models can be compared on the identical mesh. Node numbers are the MATLAB
active node numbers (`nodes(:,5)`).
"""
function mesh_from_matlab(dir::AbstractString, tag::AbstractString, inp::BasePlateInput)
    rd(name) = readdlm(joinpath(dir, "$(tag)_$(name).txt"), ',', Float64)
    nodes = rd("nodes"); el_nodes = rd("el_nodes"); el_props = rd("el_props")
    beam_nodes = rd("beam_nodes"); beam_props = rd("beam_props")
    act = Int.(nodes[:, 5])
    nn = maximum(act)
    x = zeros(nn, 3)
    for i in axes(nodes, 1)
        nodes[i, 4] > 0 && (x[act[i], :] = nodes[i, 1:3])
    end
    g = inp.geometry
    conn = Vector{Int}[]; part = Int[]; thickness = Float64[]; hole = falses(0)
    for e in axes(el_nodes, 1)
        nen = el_props[e, 2] < 20 ? 3 : 4
        push!(conn, [act[Int(el_nodes[e, k])] for k in 1:nen])
        push!(part, Int(el_props[e, 1])); push!(thickness, el_props[e, 3]); push!(hole, el_props[e, 5] < 0)
    end
    beams = [(act[Int(beam_nodes[i, 1])], act[Int(beam_nodes[i, 2])], beam_props[i, 3] > 50 ? 1 : 2) for i in axes(beam_nodes, 1)]
    bp_nodes = act[Int.(vec(rd("nods_bp")))]
    bp_area = vec(rd("nods_area"))
    up = Int.(vec(rd("upright_points"))); co = Int.(vec(rd("bp_corner")))
    bolt_nodes = act[Int.(rd("nods_spr_b")[:, 1])]
    hole_nodes = [[b[2] for b in beams if b[3] == 2 && b[1] == m] for m in bolt_nodes]
    h, b, c = centerline(g)
    frames = [Matrix(1.0I, 3, 3)]
    for w in upright_walls(g, inp.ym_ratio * b)
        len = hypot(w.p2[1] - w.p1[1], w.p2[2] - w.p1[2])
        j1 = [(w.p2[1] - w.p1[1]) / len, (w.p2[2] - w.p1[2]) / len, 0.0]; j2 = [0.0, 0.0, 1.0]
        push!(frames, hcat(j1, j2, cross(j1, j2)))
    end
    return Mesh(x, conn, part, thickness, hole, beams, bp_nodes, bp_area, act[size(nodes, 1)], bolt_nodes, hole_nodes,
                (act[up[1]], act[up[2]]), (act[co[1]], act[co[2]], act[co[3]], act[co[4]]), frames)
end

"""
    apply_imperfection(mesh, shape, amplitude) -> Mesh

Mesh with an initial geometric imperfection of the upright: the nodes above the base plate (except
the master node) are moved by the translations of `shape` (a vector of global dofs, e.g. a buckling
mode from [`buckling_modes`](@ref)), scaled so that the largest nodal movement is `amplitude`. The
sign is chosen so that the node with the largest movement moves away from the centre of the
upright cross-section (multiply `amplitude` by −1 for the other sign).
"""
function apply_imperfection(m::Mesh, shape::AbstractVector, amplitude::Real)
    x = copy(m.x)
    isbp = falses(nnodes(m)); isbp[m.bp_nodes] .= true
    nodes = [i for i in 1:nnodes(m) if !isbp[i] && i != m.master]
    d = [shape[6*(i-1)+k] for i in nodes, k in 1:3]
    mag = [norm(d[j, :]) for j in axes(d, 1)]
    jmax = argmax(mag)
    mag[jmax] > 0 || throw(ArgumentError("the imperfection shape is zero on the upright"))
    i = nodes[jmax]
    out = dot(d[jmax, 1:2], [x[i, 1] - m.x[m.master, 1], x[i, 2] - m.x[m.master, 2]])
    s = amplitude / mag[jmax] * (out < 0 ? -1 : 1)
    x[nodes, :] .+= s .* d
    return Mesh(x, m.conn, m.part, m.thickness, m.hole, m.beams, m.bp_nodes, m.bp_area, m.master, m.bolt_nodes,
                m.hole_nodes, m.rot_nodes, m.corner_nodes, m.frames)
end

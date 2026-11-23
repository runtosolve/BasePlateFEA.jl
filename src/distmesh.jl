# --------------------------------------------------------------------------------------------
# 2D DistMesh (Persson & Strang 2004), following the `distmesh.m` used by S. Ádány's column-base
# code (the FEATool variant, algorithm IALG = 2, uniform size function), with the Delaunay
# triangulation from DelaunayTriangulation.jl.
# --------------------------------------------------------------------------------------------

# unique points (tolerance `s`, by rounding as in l_deduplicate), keeping the first occurrences in order
function _dedup_points(p::Vector{NTuple{2,Float64}}, s::Float64)
    seen = Set{NTuple{2,Int}}()
    out = NTuple{2,Float64}[]
    for q in p
        key = (round(Int, q[1] / s), round(Int, q[2] / s))
        if !(key in seen)
            push!(seen, key)
            push!(out, q)
        end
    end
    return out
end

function _range_tol(p::Vector{NTuple{2,Float64}})
    xs = extrema(q[1] for q in p); ys = extrema(q[2] for q in p)
    return 1e-6 * max(xs[2] - xs[1], ys[2] - ys[1])
end

function _delaunay(p::Vector{NTuple{2,Float64}})
    tri = DelaunayTriangulation.triangulate(p; randomise = false)
    return [NTuple{3,Int}(T) for T in DelaunayTriangulation.each_solid_triangle(tri)]
end

_area2(p, t) = (p[t[2]][1] - p[t[1]][1]) * (p[t[3]][2] - p[t[1]][2]) - (p[t[2]][2] - p[t[1]][2]) * (p[t[3]][1] - p[t[1]][1])

# l_triangulate: deduplicate, triangulate, drop the triangles whose centroid is outside, orient ccw
function _triangulate(p, fd, dtrm_tol)
    p = _dedup_points(p, _range_tol(p))
    t = NTuple{3,Int}[]
    for T in _delaunay(p)
        c = ((p[T[1]][1] + p[T[2]][1] + p[T[3]][1]) / 3, (p[T[1]][2] + p[T[2]][2] + p[T[3]][2]) / 3)
        fd(c) < dtrm_tol || continue
        a = _area2(p, T) / 2
        abs(a) < 10 * eps() && continue
        push!(t, a < 0 ? (T[2], T[1], T[3]) : T)
    end
    return p, t
end

"""
    distmesh2d(fd, h0, bbox, pfix; it_max = 1000) -> (p, t)

Uniform triangular mesh of the region `fd(x) < 0` (`fd` is a signed distance function of a point
`(x, y)`) with target edge length `h0`. `bbox = ((xmin, ymin), (xmax, ymax))`; the points `pfix`
are kept as mesh nodes and come first in `p`. Returns the node coordinates and the triangles
(counter-clockwise node triplets).
"""
function distmesh2d(fd, h0::Float64, bbox, pfix::Vector{NTuple{2,Float64}}; it_max::Int = 1000)
    dp_tol = -0.001 * h0
    dtrm_tol = -0.001 * h0
    rt_tol = 0.3
    F_scale = 1.2
    F_DCF = 2.0
    dp_scale = 0.2
    dpc_tol = 0.001 * h0
    gradeps = sqrt(eps()) * h0
    IT_MIN = 20; N_RECV = 2; N_DCF = 30

    # initial points: equilateral triangle grid
    xs = bbox[1][1]:h0:bbox[2][1]
    ys = bbox[1][2]:(h0 * sqrt(3) / 2):bbox[2][2]
    p = NTuple{2,Float64}[]
    for (j, y) in enumerate(ys), x in xs
        q = (x + (iseven(j) ? h0 / 2 : 0.0), y)
        fd(q) < -dp_tol && push!(p, q)
    end
    isempty(p) && return p, NTuple{3,Int}[]

    pfix = isempty(pfix) ? pfix : _dedup_points(pfix, _range_tol(pfix))
    nfix = length(pfix)
    fixset = Set(pfix)
    p = vcat(pfix, [q for q in p if !(q in fixset)])

    t = NTuple{3,Int}[]
    edges = NTuple{2,Int}[]
    p0 = fill((Inf, Inf), length(p))
    it = 0
    while it < it_max
        it += 1
        moved = length(p0) != length(p) ? Inf : maximum(hypot(p[i][1] - p0[i][1], p[i][2] - p0[i][2]) for i in eachindex(p))
        if rt_tol * h0 < moved
            p, t = _triangulate(p, fd, dtrm_tol)
            p0 = copy(p)
            es = Set{NTuple{2,Int}}()
            for T in t, (a, b) in ((T[1], T[2]), (T[2], T[3]), (T[3], T[1]))
                push!(es, a < b ? (a, b) : (b, a))
            end
            edges = sort!(collect(es))
        end
        np = length(p)

        L = [hypot(p[a][1] - p[b][1], p[a][2] - p[b][2]) for (a, b) in edges]
        # uniform size function: every bar has the same target length
        Ltarget = F_scale * sqrt(sum(abs2, L) / length(L))

        # density control: remove the (non-fixed) end points of bars that are far too short
        if it % N_DCF == 0 && any(Ltarget .> F_DCF .* L)
            drop = Set{Int}()
            for (k, (a, b)) in enumerate(edges)
                if Ltarget > F_DCF * L[k]
                    a > nfix && push!(drop, a)
                    b > nfix && push!(drop, b)
                end
            end
            p = [p[i] for i in 1:np if !(i in drop)]
            p0 = fill((Inf, Inf), 0)
            continue
        end

        dx = zeros(np); dy = zeros(np)
        for (k, (a, b)) in enumerate(edges)
            F = max(Ltarget - L[k], 0.0)
            F == 0 && continue
            fx = F / L[k] * (p[a][1] - p[b][1]); fy = F / L[k] * (p[a][2] - p[b][2])
            dx[a] += fx; dy[a] += fy
            dx[b] -= fx; dy[b] -= fy
        end
        dx[1:nfix] .= 0; dy[1:nfix] .= 0
        dx .*= dp_scale; dy .*= dp_scale
        p = [(p[i][1] + dx[i], p[i][2] + dy[i]) for i in 1:np]

        # bring the points that left the region back to the boundary
        dist = zeros(np)
        for _ in 1:N_RECV
            for i in 1:np
                dist[i] = fd(p[i])
                if i > nfix && dist[i] > 0
                    gx = (fd((p[i][1] + gradeps, p[i][2])) - dist[i]) / gradeps
                    gy = (fd((p[i][1], p[i][2] + gradeps)) - dist[i]) / gradeps
                    g2 = gx^2 + gy^2
                    p[i] = (p[i][1] - dist[i] / g2 * gx, p[i][2] - dist[i] / g2 * gy)
                end
            end
        end

        dmax = 0.0
        for i in 1:np
            dist[i] < dp_tol && (dmax = max(dmax, hypot(dx[i], dy[i])))
        end
        ((it > IT_MIN && dmax < dpc_tol) || length(t) <= 2) && break
    end

    # l_fixmesh: final triangulation, drop the unused points
    p, t = _triangulate(p, fd, dtrm_tol)
    used = sort!(unique(vcat([collect(T) for T in t]...)))
    newid = zeros(Int, length(p))
    newid[used] = 1:length(used)
    return p[used], [(newid[T[1]], newid[T[2]], newid[T[3]]) for T in t]
end

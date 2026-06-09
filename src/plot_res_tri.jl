function plot_res_tri(el_nodes, el_props, out_grid, FE_dat, out_res, scale, layer, comp, gridon;
                      title_str = "Component $comp",
                      colormap  = :jet,
                      azimuth   = deg2rad(210),
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

    # ── figure ────────────────────────────────────────────────────────────────
    fig = Figure(size = (1100, 800))
    ax  = Axis3(fig[1, 1];
                aspect    = :data,
                xlabel    = "X (in)",
                ylabel    = "Y (in)",
                zlabel    = "Z (in)",
                title     = title_str,
                azimuth   = azimuth,
                elevation = elevation)

    # filled colour-mapped mesh with smooth per-vertex interpolation
    msh = mesh!(ax, geom;
                color    = out_comp,
                colormap = colormap,
                shading  = NoShading)

    # optional grey wireframe edge overlay
    if gridon
        wireframe!(ax, geom;
                   color     = (:gray, 0.35),
                   linewidth = 0.5)
    end

    Colorbar(fig[1, 2], msh; width = 20, label = title_str)

    return fig
end

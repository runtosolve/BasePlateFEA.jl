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

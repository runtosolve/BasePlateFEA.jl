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
    LV = zeros(totaldof)
    for k in 1:i0
        LV[rows[k]] += vals[k]
    end
    return LV
end

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

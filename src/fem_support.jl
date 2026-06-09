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

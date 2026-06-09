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

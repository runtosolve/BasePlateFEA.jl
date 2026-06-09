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

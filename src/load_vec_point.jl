function load_vec_point(nodes, nodenr, PLATE_table, LOAD_table, dpn, tolers)
    tol = tolers.zero_dist
    totaldof = nodenr * dpn
    LV = zeros(totaldof)

    for i in 1:length(LOAD_table)
        LOAD_table[i][1] != "point" && continue
        ipl = LOAD_table[i][2]
        L1p = [Float64(LOAD_table[i][3]), Float64(LOAD_table[i][4]), 0.0]
        px = Float64(LOAD_table[i][5]); py = Float64(LOAD_table[i][6]); pz = Float64(LOAD_table[i][7])
        mx = Float64(LOAD_table[i][8]); my = Float64(LOAD_table[i][9]); mz = Float64(LOAD_table[i][10])
        P1 = PLATE_table[ipl][1:3]; P2 = PLATE_table[ipl][4:6]
        P3 = PLATE_table[ipl][7:9]; P4 = PLATE_table[ipl][10:12]
        _, T3 = plate_transform(P1, P2, P3, P4)
        L1 = T3 * L1p + P1
        pvec = [px, py, pz, mx, my, mz]
        ltype_flag = LOAD_table[i][11]
        if isa(ltype_flag, String) && startswith(ltype_flag, "s")
            T6 = rotate3d(T3, 1, 6)
            pvec = T6 * pvec
        end
        dif = nodes[:, 1:3] .- L1'
        ind = findall(j -> abs(dif[j,1]) < tol && abs(dif[j,2]) < tol && abs(dif[j,3]) < tol, 1:size(nodes,1))
        isempty(ind) && continue
        active_nod = nodes[ind[1], 5]
        LV[(active_nod-1)*dpn+1 : (active_nod-1)*dpn+dpn] .+= pvec[1:dpn]
    end
    return LV
end

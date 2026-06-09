function collect_anchor(PLATE_table, SUPPORT_table, CONDSUP_table, tol)
    npl = length(PLATE_table)
    p_anchor = [zeros(0, 3) for _ in 1:npl]

    # anchor points from plate geometry intersections
    for ipl in 1:npl
        Q = zeros(4, 3)
        Q[1, :] = PLATE_table[ipl][1:3]
        Q[2, :] = PLATE_table[ipl][4:6]
        Q[3, :] = PLATE_table[ipl][7:9]
        Q[4, :] = PLATE_table[ipl][10:12]
        for iside in 1:4
            i1 = iside; i2 = iside == 4 ? 1 : iside + 1
            Q1 = Q[i1, :]; Q2 = Q[i2, :]
            for ipl2 in 1:npl
                ipl2 == ipl && continue
                for icorner in 1:4
                    idx = (icorner-1)*3+1 : (icorner-1)*3+3
                    P = PLATE_table[ipl2][idx]
                    if abs(Q1[2]-Q2[2]) < tol && abs(Q1[3]-Q2[3]) < tol && P[1] > min(Q1[1],Q2[1])+tol && P[1] < max(Q1[1],Q2[1])-tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], [P[1] Q1[2] Q1[3]])
                    end
                    if abs(Q1[1]-Q2[1]) < tol && abs(Q1[3]-Q2[3]) < tol && P[2] > min(Q1[2],Q2[2])+tol && P[2] < max(Q1[2],Q2[2])-tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], [Q1[1] P[2] Q1[3]])
                    end
                    if abs(Q1[1]-Q2[1]) < tol && abs(Q1[2]-Q2[2]) < tol && P[3] > min(Q1[3],Q2[3])+tol && P[3] < max(Q1[3],Q2[3])-tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], [Q1[1] Q1[2] P[3]])
                    end
                end
            end
        end
    end

    # anchor points from ordinary supports
    for ipl in 1:npl
        Q = zeros(4, 3)
        Q[1, :] = PLATE_table[ipl][1:3]
        Q[2, :] = PLATE_table[ipl][4:6]
        Q[3, :] = PLATE_table[ipl][7:9]
        Q[4, :] = PLATE_table[ipl][10:12]
        for iside in 1:4
            i1 = iside; i2 = iside == 4 ? 1 : iside + 1
            Q1 = Q[i1, :]; Q2 = Q[i2, :]
            for isup in 1:length(SUPPORT_table)
                stype = SUPPORT_table[isup][1]
                if stype in ("point", "p-line", "p-surf", "line", "surf")
                    jpl = SUPPORT_table[isup][2]
                    P1j = PLATE_table[jpl][1:3]; P2j = PLATE_table[jpl][4:6]
                    P3j = PLATE_table[jpl][7:9]; P4j = PLATE_table[jpl][10:12]
                    _, T = plate_transform(P1j, P2j, P3j, P4j)
                    if stype in ("point",)
                        npoi = 1; indpoi = [(3,4)]
                    else
                        npoi = 2; indpoi = [(3,4), (5,6)]
                    end
                    for ipoi in 1:npoi
                        ii1, ii2 = indpoi[ipoi]
                        Pp = [Float64(SUPPORT_table[isup][ii1]), Float64(SUPPORT_table[isup][ii2]), 0.0]
                        P = T * Pp + P1j
                        if abs(Q1[2]-Q2[2]) < tol && abs(Q1[3]-Q2[3]) < tol && P[1] > min(Q1[1],Q2[1])+tol && P[1] < max(Q1[1],Q2[1])-tol
                            p_anchor[ipl] = vcat(p_anchor[ipl], [P[1] Q1[2] Q1[3]])
                        end
                        if abs(Q1[1]-Q2[1]) < tol && abs(Q1[3]-Q2[3]) < tol && P[2] > min(Q1[2],Q2[2])+tol && P[2] < max(Q1[2],Q2[2])-tol
                            p_anchor[ipl] = vcat(p_anchor[ipl], [Q1[1] P[2] Q1[3]])
                        end
                        if abs(Q1[1]-Q2[1]) < tol && abs(Q1[2]-Q2[2]) < tol && P[3] > min(Q1[3],Q2[3])+tol && P[3] < max(Q1[3],Q2[3])-tol
                            p_anchor[ipl] = vcat(p_anchor[ipl], [Q1[1] Q1[2] P[3]])
                        end
                    end
                end
            end
        end
    end

    # anchor points from conditional supports
    for ipl in 1:npl
        Q = zeros(4, 3)
        Q[1, :] = PLATE_table[ipl][1:3]
        Q[2, :] = PLATE_table[ipl][4:6]
        Q[3, :] = PLATE_table[ipl][7:9]
        Q[4, :] = PLATE_table[ipl][10:12]
        for iside in 1:4
            i1 = iside; i2 = iside == 4 ? 1 : iside + 1
            Q1 = Q[i1, :]; Q2 = Q[i2, :]
            for isup in 1:length(CONDSUP_table)
                stype = CONDSUP_table[isup][1]
                jpl = CONDSUP_table[isup][2]
                P1j = PLATE_table[jpl][1:3]; P2j = PLATE_table[jpl][4:6]
                P3j = PLATE_table[jpl][7:9]; P4j = PLATE_table[jpl][10:12]
                _, T = plate_transform(P1j, P2j, P3j, P4j)
                pts_to_check = Vector{Vector{Float64}}()
                if stype == "point"
                    push!(pts_to_check, [Float64(CONDSUP_table[isup][3]), Float64(CONDSUP_table[isup][4]), 0.0])
                elseif stype == "surf"
                    push!(pts_to_check, [Float64(CONDSUP_table[isup][3]), Float64(CONDSUP_table[isup][4]), 0.0])
                    push!(pts_to_check, [Float64(CONDSUP_table[isup][5]), Float64(CONDSUP_table[isup][6]), 0.0])
                end
                for Pp in pts_to_check
                    P = T * Pp + P1j
                    if abs(Q1[2]-Q2[2]) < tol && abs(Q1[3]-Q2[3]) < tol && P[1] > min(Q1[1],Q2[1])+tol && P[1] < max(Q1[1],Q2[1])-tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], [P[1] Q1[2] Q1[3]])
                    end
                    if abs(Q1[1]-Q2[1]) < tol && abs(Q1[3]-Q2[3]) < tol && P[2] > min(Q1[2],Q2[2])+tol && P[2] < max(Q1[2],Q2[2])-tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], [Q1[1] P[2] Q1[3]])
                    end
                    if abs(Q1[1]-Q2[1]) < tol && abs(Q1[2]-Q2[2]) < tol && P[3] > min(Q1[3],Q2[3])+tol && P[3] < max(Q1[3],Q2[3])-tol
                        p_anchor[ipl] = vcat(p_anchor[ipl], [Q1[1] Q1[2] P[3]])
                    end
                end
            end
        end
    end

    return p_anchor
end

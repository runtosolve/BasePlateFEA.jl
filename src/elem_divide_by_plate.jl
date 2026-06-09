function elem_divide_by_plate(pl_elems, pl_nodes, PLATE_table, HOLE_table, tol)
    for ipl in 1:length(pl_elems)
        HOLE_TABLE_ipl = HOLE_TABLE_entry(HOLE_table, ipl)
        isnothing(HOLE_TABLE_ipl) && continue

        C1 = PLATE_table[ipl][1:3]; C2 = PLATE_table[ipl][4:6]
        C3 = PLATE_table[ipl][7:9]; C4 = PLATE_table[ipl][10:12]
        _, T = plate_transform(C1, C2, C3, C4)

        el_nodes = pl_elems[ipl]
        nelems = size(el_nodes, 1)
        eln2 = zeros(Int, 0, 9)
        elim_nodes = Int[]

        for ihole in 1:length(HOLE_TABLE_ipl.holes)
            npoly = size(HOLE_TABLE_ipl.holes[ihole], 1) - 1
            for iel in 1:nelems
                found = false; ipoly = 0
                while !found && ipoly < npoly
                    ipoly += 1
                    nsec = 0
                    P1h = HOLE_TABLE_ipl.holes[ihole][ipoly,   1:2]
                    P2h = HOLE_TABLE_ipl.holes[ihole][ipoly+1, 1:2]
                    el_sec_points = zeros(4, 4)
                    for iside in 1:4
                        i1 = iside; i2 = iside == 4 ? 1 : iside + 1
                        Q1g = pl_nodes[ipl][el_nodes[iel, i1], 1:3]
                        Q2g = pl_nodes[ipl][el_nodes[iel, i2], 1:3]
                        Q1 = (T' * (Q1g - C1))[1:2]
                        Q2 = (T' * (Q2g - C1))[1:2]
                        A = [P2h[1]-P1h[1]  Q1[1]-Q2[1];
                             P2h[2]-P1h[2]  Q1[2]-Q2[2]]
                        b_vec = [Q1[1]-P1h[1]; Q1[2]-P1h[2]]
                        if abs(det(A)) > tol
                            x = A \ b_vec
                            if x[1] >= -tol && x[1] <= 1+tol && x[2] >= -tol && x[2] <= 1+tol
                                if (x[1] >= tol && x[1] <= 1-tol) || (x[2] >= tol && x[2] <= 1-tol)
                                    # mesh error - just note it
                                else
                                    nsec += 1
                                    el_sec_points[nsec, :] = [ipoly, iel, iside, round(x[2])]
                                end
                            end
                        end
                    end
                    if nsec == 4
                        new_row = zeros(Int, 1, 9)
                        if el_sec_points[1, 4] == 1
                            new_row[1, 1:3] = el_nodes[iel, [3, 4, 2]]
                            el_nodes[iel, 1:3] = el_nodes[iel, [1, 2, 4]]
                        else
                            new_row[1, 1:3] = el_nodes[iel, [4, 1, 3]]
                            el_nodes[iel, 1:3] = el_nodes[iel, [2, 3, 1]]
                        end
                        eln2 = vcat(eln2, new_row)
                        push!(elim_nodes, iel)
                        found = true
                    end
                end
            end
        end

        for iel in elim_nodes
            el_nodes[iel, 4] = 0
        end
        pl_elems[ipl] = vcat(el_nodes, eln2)
    end
    return pl_elems
end

function HOLE_TABLE_entry(HOLE_table, ipl)
    if ipl > length(HOLE_table) || isnothing(HOLE_table[ipl])
        return nothing
    end
    entry = HOLE_table[ipl]
    if isnothing(entry) || (isa(entry, HoleData) && isempty(entry.holes))
        return nothing
    end
    return entry
end

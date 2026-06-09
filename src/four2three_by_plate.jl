function four2three_by_plate(pl_elems, rat_q2t)
    for ipl in 1:length(pl_elems)
        el_nodes = pl_elems[ipl]
        nelems = size(el_nodes, 1)
        eln2 = zeros(Int, 0, 9)
        for iel in 1:nelems
            if el_nodes[iel, 4] > 0 && rand() > rat_q2t[ipl]
                new_row = zeros(Int, 1, 9)
                if rand() > (1 - rat_q2t[ipl]) / 2 + rat_q2t[ipl]
                    new_row[1, 1:3] = el_nodes[iel, [3, 4, 2]]
                    el_nodes[iel, 1:3] = el_nodes[iel, [1, 2, 4]]
                else
                    new_row[1, 1:3] = el_nodes[iel, [4, 1, 3]]
                    el_nodes[iel, 1:3] = el_nodes[iel, [2, 3, 1]]
                end
                el_nodes[iel, 4] = 0
                eln2 = vcat(eln2, new_row)
            end
        end
        pl_elems[ipl] = vcat(el_nodes, eln2)
    end
    return pl_elems
end

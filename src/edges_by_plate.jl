function edges_by_plate(pl_elems)
    pl_el_edges = Vector{Matrix{Int}}()
    pl_edge_nodes = Vector{Matrix{Int}}()

    for ipl in 1:length(pl_elems)
        el_nodes = pl_elems[ipl]
        nelems = size(el_nodes, 1)
        edge_tot = 0
        edge_nodes = zeros(Int, 0, 2)
        el_edges = zeros(Int, nelems, 4)

        for iel in 1:nelems
            nedg = el_nodes[iel, 4] == 0 ? 3 : 4
            for iedg in 1:nedg
                node1 = el_nodes[iel, iedg]
                node2 = iedg < nedg ? el_nodes[iel, iedg+1] : el_nodes[iel, 1]
                ind = findall(i -> (edge_nodes[i,1]==node1 && edge_nodes[i,2]==node2) ||
                                   (edge_nodes[i,1]==node2 && edge_nodes[i,2]==node1),
                              1:size(edge_nodes,1))
                if !isempty(ind)
                    el_edges[iel, iedg] = ind[1]
                else
                    edge_tot += 1
                    el_edges[iel, iedg] = edge_tot
                    edge_nodes = vcat(edge_nodes, [node1 node2])
                end
            end
        end
        push!(pl_el_edges, el_edges)
        push!(pl_edge_nodes, edge_nodes)
    end
    return pl_el_edges, pl_edge_nodes
end

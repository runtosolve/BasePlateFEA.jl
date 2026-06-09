function ct_3node_g2e(elems_row::AbstractVector, nodes::AbstractMatrix)
    i1 = elems_row[1]; i2 = elems_row[2]; i3 = elems_row[3]
    P1g = nodes[i1, 1:3]; P2g = nodes[i2, 1:3]; P3g = nodes[i3, 1:3]
    norm_vec = cross(P2g - P1g, P3g - P1g)
    norm_vec = norm_vec / norm(norm_vec)
    j3 = norm_vec
    j1 = (P2g - P1g) / norm(P2g - P1g)
    j2 = cross(j3, j1)
    T = hcat(j1, j2, j3)   # 3×3, columns are local axes
    P1e = T' * (P1g - P1g)
    P2e = T' * (P2g - P1g)
    P3e = T' * (P3g - P1g)
    return T, P1e, P2e, P3e
end

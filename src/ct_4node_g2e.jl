function ct_4node_g2e(elems_row::AbstractVector, nodes::AbstractMatrix)
    i1 = elems_row[1]; i2 = elems_row[2]; i3 = elems_row[3]; i4 = elems_row[4]
    P1g = nodes[i1, 1:3]; P2g = nodes[i2, 1:3]
    P3g = nodes[i3, 1:3]; P4g = nodes[i4, 1:3]
    P12 = (P1g + P2g) / 2;  P34 = (P3g + P4g) / 2
    P23 = (P2g + P3g) / 2;  P41 = (P4g + P1g) / 2
    norm_vec = cross(P23 - P41, P34 - P12)
    norm_vec = norm_vec / norm(norm_vec)
    j3 = norm_vec
    P0 = P12 + (P34 - P12) / 2

    P1 = (P1g - P0) - dot(P1g - P0, norm_vec) * norm_vec
    P2 = (P2g - P0) - dot(P2g - P0, norm_vec) * norm_vec
    P3 = (P3g - P0) - dot(P3g - P0, norm_vec) * norm_vec
    P4 = (P4g - P0) - dot(P4g - P0, norm_vec) * norm_vec

    P12 = (P1 + P2) / 2;  P34 = (P3 + P4) / 2
    P23 = (P2 + P3) / 2;  P41 = (P4 + P1) / 2
    j1 = (P23 - P41) / norm(P23 - P41)
    j2 = cross(j3, j1)
    T = hcat(j1, j2, j3)
    P1e = T' * P1
    P2e = T' * P2
    P3e = T' * P3
    P4e = T' * P4
    return T, P1e, P2e, P3e, P4e
end

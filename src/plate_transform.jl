function plate_transform(P1::AbstractVector, P2::AbstractVector,
                         P3::AbstractVector, P4::AbstractVector)
    P12 = (P1 + P2) / 2;  P34 = (P3 + P4) / 2
    P23 = (P2 + P3) / 2;  P41 = (P4 + P1) / 2
    norm_vec = cross(P23 - P41, P34 - P12)
    norm_vec = norm_vec / norm(norm_vec)
    j3 = norm_vec
    P0 = P12 + (P34 - P12) / 2

    P1p = (P1 - P0) - dot(P1 - P0, norm_vec) * norm_vec
    P2p = (P2 - P0) - dot(P2 - P0, norm_vec) * norm_vec
    P3p = (P3 - P0) - dot(P3 - P0, norm_vec) * norm_vec
    P4p = (P4 - P0) - dot(P4 - P0, norm_vec) * norm_vec

    P12 = (P1p + P2p) / 2;  P34 = (P3p + P4p) / 2
    P23 = (P2p + P3p) / 2;  P41 = (P4p + P1p) / 2
    i23_41 = (P23 - P41) / norm(P23 - P41)
    j1 = i23_41
    j2 = cross(j3, j1)

    T = hcat(j1, j2, j3)   # 3×3, columns are local axes
    return P0, T
end

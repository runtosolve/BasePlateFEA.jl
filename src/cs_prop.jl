function cs_prop(cs_node::AbstractMatrix, cs_elem::AbstractMatrix)
    np = size(cs_elem, 1)
    Y1 = zeros(np); Z1 = zeros(np)
    Y2 = zeros(np); Z2 = zeros(np)
    b_width = zeros(np); alph = zeros(np); t = zeros(np)

    for i in 1:np
        nodei = Int(cs_elem[i, 2])
        nodej = Int(cs_elem[i, 3])
        yi = cs_node[nodei, 2]; zi = cs_node[nodei, 3]
        yj = cs_node[nodej, 2]; zj = cs_node[nodej, 3]
        Y1[i] = yi; Z1[i] = zi
        Y2[i] = yj; Z2[i] = zj
        dy = yj - yi; dz = zj - zi
        width = sqrt(dy^2 + dz^2)
        b_width[i] = width
        alph[i] = atan(dz, dy)
        t[i] = cs_elem[i, 4]
    end
    Ym = (Y1 + Y2) / 2
    Zm = (Z1 + Z2) / 2

    dA = b_width .* t
    A = sum(dA)

    Sy = sum(Zm .* dA)
    Sz = sum(Ym .* dA)
    yCG = Sz / A
    zCG = Sy / A

    Iy = sum(t .* ((sin.(alph).^2 .* b_width.^2 + cos.(alph).^2 .* t.^2) / 12 + Zm.^2) .* b_width) - A * zCG^2
    Iz = sum(t .* ((cos.(alph).^2 .* b_width.^2 + sin.(alph).^2 .* t.^2) / 12 + Ym.^2) .* b_width) - A * yCG^2

    Iyr = sum((1/3 .* (Z1.^2 + Z2.^2 + Z1.*Z2)) .* dA) - A * zCG^2
    Izr = sum((1/3 .* (Y1.^2 + Y2.^2 + Y1.*Y2)) .* dA) - A * yCG^2

    return A, yCG, zCG, Iyr, Iy, Izr, Iz
end

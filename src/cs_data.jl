function cs_data(h, b, c, tb, th, tc)
    cs_node, cs_elem = cs_def(h, b, c, tb, th, tc)
    A, yCG, zCG, Iyr, Iy, Izr, Iz = cs_prop(cs_node, cs_elem)
    return A, yCG, zCG, Iyr, Iy, Izr, Iz
end

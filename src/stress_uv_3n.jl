function stress_uv_3n(d_elem, P1, P2, P3, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp, ifgauss)
    xe1=P1[1]; ye1=P1[2]
    xe2=P2[1]; ye2=P2[2]
    xe3=P3[1]; ye3=P3[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)

    D = xe1*ye2 - xe2*ye1 - xe1*ye3 + xe3*ye1 + xe2*ye3 - xe3*ye2

    DB = [
        (E11*t*(ye2-ye3))/D  -(E12*t*(xe2-xe3))/D  -(E11*t*(ye1-ye3))/D  (E12*t*(xe1-xe3))/D  (E11*t*(ye1-ye2))/D  -(E12*t*(xe1-xe2))/D;
        (E21*t*(ye2-ye3))/D  -(E22*t*(xe2-xe3))/D  -(E21*t*(ye1-ye3))/D  (E22*t*(xe1-xe3))/D  (E21*t*(ye1-ye2))/D  -(E22*t*(xe1-xe2))/D;
        -(G*t*(xe2-xe3))/D    (G*t*(ye2-ye3))/D     (G*t*(xe1-xe3))/D   -(G*t*(ye1-ye3))/D   -(G*t*(xe1-xe2))/D    (G*t*(ye1-ye2))/D
    ]

    npoint = ifgauss ? nGp : 3
    str_uv = zeros(3, max(npoint, 1))
    str_uv[:, 1] = DB * d_elem / t
    for i in 2:npoint
        str_uv[:, i] = str_uv[:, 1]
    end
    return str_uv
end

function ke_uv_3n_condens_from_8to6dof(P1, P2, P3, t, Ex, Ey, nuxy, nuyx, G, mat_sym)
    xe1=P1[1]; ye1=P1[2]
    xe2=P2[1]; ye2=P2[2]
    xe3=P3[1]; ye3=P3[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)

    D = xe1*ye2 - xe2*ye1 - xe1*ye3 + xe3*ye1 + xe2*ye3 - xe3*ye2

    kuv = [
        (G*t*xe2^2-2*G*t*xe2*xe3+G*t*xe3^2+E11*t*ye2^2-2*E11*t*ye2*ye3+E11*t*ye3^2)/(2*D)  -(t*(xe2-xe3)*(ye2-ye3)*(E12+G))/(2*D)  -(E11*t*ye3^2+G*t*xe3^2+G*t*xe1*xe2-G*t*xe1*xe3-G*t*xe2*xe3+E11*t*ye1*ye2-E11*t*ye1*ye3-E11*t*ye2*ye3)/(2*D)  (E12*t*xe1*ye2-E12*t*xe1*ye3-E12*t*xe3*ye2+E12*t*xe3*ye3+G*t*xe2*ye1-G*t*xe3*ye1-G*t*xe2*ye3+G*t*xe3*ye3)/(2*D)  -(E11*t*ye2^2+G*t*xe2^2-G*t*xe1*xe2+G*t*xe1*xe3-G*t*xe2*xe3-E11*t*ye1*ye2+E11*t*ye1*ye3-E11*t*ye2*ye3)/(2*D)  -(E12*t*xe1*ye2-E12*t*xe1*ye3-E12*t*xe2*ye2+E12*t*xe2*ye3+G*t*xe2*ye1-G*t*xe2*ye2-G*t*xe3*ye1+G*t*xe3*ye2)/(2*D)  0.0  0.0
        -(t*(xe2-xe3)*(ye2-ye3)*(E21+G))/(2*D)  (E22*t*xe2^2-2*E22*t*xe2*xe3+E22*t*xe3^2+G*t*ye2^2-2*G*t*ye2*ye3+G*t*ye3^2)/(2*D)  (E21*t*xe2*ye1-E21*t*xe3*ye1-E21*t*xe2*ye3+E21*t*xe3*ye3+G*t*xe1*ye2-G*t*xe1*ye3-G*t*xe3*ye2+G*t*xe3*ye3)/(2*D)  -(E22*t*xe3^2+G*t*ye3^2+E22*t*xe1*xe2-E22*t*xe1*xe3-E22*t*xe2*xe3+G*t*ye1*ye2-G*t*ye1*ye3-G*t*ye2*ye3)/(2*D)  -(E21*t*xe2*ye1-E21*t*xe2*ye2-E21*t*xe3*ye1+E21*t*xe3*ye2+G*t*xe1*ye2-G*t*xe1*ye3-G*t*xe2*ye2+G*t*xe2*ye3)/(2*D)  -(E22*t*xe2^2+G*t*ye2^2-E22*t*xe1*xe2+E22*t*xe1*xe3-E22*t*xe2*xe3-G*t*ye1*ye2+G*t*ye1*ye3-G*t*ye2*ye3)/(2*D)  0.0  0.0
        -(E11*t*ye3^2+G*t*xe3^2+G*t*xe1*xe2-G*t*xe1*xe3-G*t*xe2*xe3+E11*t*ye1*ye2-E11*t*ye1*ye3-E11*t*ye2*ye3)/(2*D)  (E12*t*xe2*ye1-E12*t*xe3*ye1-E12*t*xe2*ye3+E12*t*xe3*ye3+G*t*xe1*ye2-G*t*xe1*ye3-G*t*xe3*ye2+G*t*xe3*ye3)/(2*D)  (t*(G*xe1^2-2*G*xe1*xe3+G*xe3^2+E11*ye1^2-2*E11*ye1*ye3+E11*ye3^2))/(2*D)  -(t*(xe1-xe3)*(ye1-ye3)*(E12+G))/(2*D)  -(t*(E11*ye1^2+G*xe1^2-G*xe1*xe2-G*xe1*xe3+G*xe2*xe3-E11*ye1*ye2-E11*ye1*ye3+E11*ye2*ye3))/(2*D)  (t*(E12*xe1*ye1-E12*xe2*ye1-E12*xe1*ye3+E12*xe2*ye3+G*xe1*ye1-G*xe1*ye2-G*xe3*ye1+G*xe3*ye2))/(2*D)  0.0  0.0
        (E21*t*xe1*ye2-E21*t*xe1*ye3-E21*t*xe3*ye2+E21*t*xe3*ye3+G*t*xe2*ye1-G*t*xe3*ye1-G*t*xe2*ye3+G*t*xe3*ye3)/(2*D)  -(E22*t*xe3^2+G*t*ye3^2+E22*t*xe1*xe2-E22*t*xe1*xe3-E22*t*xe2*xe3+G*t*ye1*ye2-G*t*ye1*ye3-G*t*ye2*ye3)/(2*D)  -(t*(xe1-xe3)*(ye1-ye3)*(E21+G))/(2*D)  (t*(E22*xe1^2-2*E22*xe1*xe3+E22*xe3^2+G*ye1^2-2*G*ye1*ye3+G*ye3^2))/(2*D)  (t*(E21*xe1*ye1-E21*xe1*ye2-E21*xe3*ye1+E21*xe3*ye2+G*xe1*ye1-G*xe2*ye1-G*xe1*ye3+G*xe2*ye3))/(2*D)  -(t*(E22*xe1^2+G*ye1^2-E22*xe1*xe2-E22*xe1*xe3+E22*xe2*xe3-G*ye1*ye2-G*ye1*ye3+G*ye2*ye3))/(2*D)  0.0  0.0
        -(E11*t*ye2^2+G*t*xe2^2-G*t*xe1*xe2+G*t*xe1*xe3-G*t*xe2*xe3-E11*t*ye1*ye2+E11*t*ye1*ye3-E11*t*ye2*ye3)/(2*D)  -(E12*t*xe2*ye1-E12*t*xe2*ye2-E12*t*xe3*ye1+E12*t*xe3*ye2+G*t*xe1*ye2-G*t*xe1*ye3-G*t*xe2*ye2+G*t*xe2*ye3)/(2*D)  -(t*(E11*ye1^2+G*xe1^2-G*xe1*xe2-G*xe1*xe3+G*xe2*xe3-E11*ye1*ye2-E11*ye1*ye3+E11*ye2*ye3))/(2*D)  (t*(E12*xe1*ye1-E12*xe1*ye2-E12*xe3*ye1+E12*xe3*ye2+G*xe1*ye1-G*xe2*ye1-G*xe1*ye3+G*xe2*ye3))/(2*D)  (t*(G*xe1^2-2*G*xe1*xe2+G*xe2^2+E11*ye1^2-2*E11*ye1*ye2+E11*ye2^2))/(2*D)  -(t*(xe1-xe2)*(ye1-ye2)*(E12+G))/(2*D)  0.0  0.0
        -(E21*t*xe1*ye2-E21*t*xe1*ye3-E21*t*xe2*ye2+E21*t*xe2*ye3+G*t*xe2*ye1-G*t*xe2*ye2-G*t*xe3*ye1+G*t*xe3*ye2)/(2*D)  -(E22*t*xe2^2+G*t*ye2^2-E22*t*xe1*xe2+E22*t*xe1*xe3-E22*t*xe2*xe3-G*t*ye1*ye2+G*t*ye1*ye3-G*t*ye2*ye3)/(2*D)  (t*(E21*xe1*ye1-E21*xe2*ye1-E21*xe1*ye3+E21*xe2*ye3+G*xe1*ye1-G*xe1*ye2-G*xe3*ye1+G*xe3*ye2))/(2*D)  -(t*(E22*xe1^2+G*ye1^2-E22*xe1*xe2-E22*xe1*xe3+E22*xe2*xe3-G*ye1*ye2-G*ye1*ye3+G*ye2*ye3))/(2*D)  -(t*(xe1-xe2)*(ye1-ye2)*(E21+G))/(2*D)  (t*(E22*xe1^2-2*E22*xe1*xe2+E22*xe2^2+G*ye1^2-2*G*ye1*ye2+G*ye2^2))/(2*D)  0.0  0.0
        0.0  0.0  0.0  0.0  0.0  0.0  (81*t*(G*xe1^2-G*xe1*xe2-G*xe1*xe3+G*xe2^2-G*xe2*xe3+G*xe3^2+E11*ye1^2-E11*ye1*ye2-E11*ye1*ye3+E11*ye2^2-E11*ye2*ye3+E11*ye3^2))/(20*D)  (81*t*(E12+G)*(xe1*ye2-2*xe1*ye1+xe2*ye1+xe1*ye3-2*xe2*ye2+xe3*ye1+xe2*ye3+xe3*ye2-2*xe3*ye3))/(40*D)
        0.0  0.0  0.0  0.0  0.0  0.0  (81*t*(E21+G)*(xe1*ye2-2*xe1*ye1+xe2*ye1+xe1*ye3-2*xe2*ye2+xe3*ye1+xe2*ye3+xe3*ye2-2*xe3*ye3))/(40*D)  (81*t*(E22*xe1^2-E22*xe1*xe2-E22*xe1*xe3+E22*xe2^2-E22*xe2*xe3+E22*xe3^2+G*ye1^2-G*ye1*ye2-G*ye1*ye3+G*ye2^2-G*ye2*ye3+G*ye3^2))/(20*D)
    ]

    # static condensation
    inda = 1:6; indi = 7:8
    ke_uv = kuv[inda, inda] - kuv[inda, indi] * (kuv[indi, indi] \ kuv[indi, inda])
    return ke_uv
end

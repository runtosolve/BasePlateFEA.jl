function ke_wt_3n_condens_from_12to9dof_MIN(P1, P2, P3, t, Ex, Ey, nuxy, nuyx, G, mat_sym)
    xe1=P1[1]; ye1=P1[2]
    xe2=P2[1]; ye2=P2[2]
    xe3=P3[1]; ye3=P3[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)

    D = xe1*ye2 - xe2*ye1 - xe1*ye3 + xe3*ye1 + xe2*ye3 - xe3*ye2

    # shear part (12×12)
    kwts = [
        (5*G*t*xe2^2-10*G*t*xe2*xe3+5*G*t*xe3^2+5*G*t*ye2^2-10*G*t*ye2*ye3+5*G*t*ye3^2)/(12*D)  (5*G*t*(xe2-xe3))/36  (5*G*t*(ye2-ye3))/36  -(5*G*t*xe3^2+5*G*t*ye3^2+5*G*t*xe1*xe2-5*G*t*xe1*xe3-5*G*t*xe2*xe3+5*G*t*ye1*ye2-5*G*t*ye1*ye3-5*G*t*ye2*ye3)/(12*D)  (5*G*t*(xe2-xe3))/36  (5*G*t*(ye2-ye3))/36  -(5*G*t*xe2^2+5*G*t*ye2^2-5*G*t*xe1*xe2+5*G*t*xe1*xe3-5*G*t*xe2*xe3-5*G*t*ye1*ye2+5*G*t*ye1*ye3-5*G*t*ye2*ye3)/(12*D)  (5*G*t*(xe2-xe3))/36  (5*G*t*(ye2-ye3))/36  -(5*G*t*(xe1*xe2-xe1*xe3+xe2*xe3+ye1*ye2-ye1*ye3+ye2*ye3-xe2^2-ye2^2))/(9*D)  -(5*G*t*(xe2^2-2*xe2*xe3+xe3^2+ye2^2-2*ye2*ye3+ye3^2))/(9*D)  (5*G*t*(xe1*xe2-xe1*xe3-xe2*xe3+ye1*ye2-ye1*ye3-ye2*ye3+xe3^2+ye3^2))/(9*D)
        (5*G*t*(xe2-xe3))/36  (5*G*t*D)/72  0.0  -(5*G*t*(xe1-xe3))/36  (5*G*t*D)/144  0.0  (5*G*t*(xe1-xe2))/36  (5*G*t*D)/144  0.0  (5*G*t*(xe2-2*xe1+xe3))/36  -(5*G*t*(xe2-xe3))/36  -(5*G*t*(xe2-2*xe1+xe3))/36
        (5*G*t*(ye2-ye3))/36  0.0  (5*G*t*D)/72  -(5*G*t*(ye1-ye3))/36  0.0  (5*G*t*D)/144  (5*G*t*(ye1-ye2))/36  0.0  (5*G*t*D)/144  (5*G*t*(ye2-2*ye1+ye3))/36  -(5*G*t*(ye2-ye3))/36  -(5*G*t*(ye2-2*ye1+ye3))/36
        -(5*G*t*xe3^2+5*G*t*ye3^2+5*G*t*xe1*xe2-5*G*t*xe1*xe3-5*G*t*xe2*xe3+5*G*t*ye1*ye2-5*G*t*ye1*ye3-5*G*t*ye2*ye3)/(12*D)  -(5*G*t*(xe1-xe3))/36  -(5*G*t*(ye1-ye3))/36  (5*G*t*(xe1^2-2*xe1*xe3+xe3^2+ye1^2-2*ye1*ye3+ye3^2))/(12*D)  -(5*G*t*(xe1-xe3))/36  -(5*G*t*(ye1-ye3))/36  (5*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(12*D)  -(5*G*t*(xe1-xe3))/36  -(5*G*t*(ye1-ye3))/36  -(5*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(9*D)  (5*G*t*(xe1*xe2-xe1*xe3-xe2*xe3+ye1*ye2-ye1*ye3-ye2*ye3+xe3^2+ye3^2))/(9*D)  -(5*G*t*(xe1^2-2*xe1*xe3+xe3^2+ye1^2-2*ye1*ye3+ye3^2))/(9*D)
        (5*G*t*(xe2-xe3))/36  (5*G*t*D)/144  0.0  -(5*G*t*(xe1-xe3))/36  (5*G*t*D)/72  0.0  (5*G*t*(xe1-xe2))/36  (5*G*t*D)/144  0.0  -(5*G*t*(xe1-2*xe2+xe3))/36  (5*G*t*(xe1-2*xe2+xe3))/36  (5*G*t*(xe1-xe3))/36
        (5*G*t*(ye2-ye3))/36  0.0  (5*G*t*D)/144  -(5*G*t*(ye1-ye3))/36  0.0  (5*G*t*D)/72  (5*G*t*(ye1-ye2))/36  0.0  (5*G*t*D)/144  -(5*G*t*(ye1-2*ye2+ye3))/36  (5*G*t*(ye1-2*ye2+ye3))/36  (5*G*t*(ye1-ye3))/36
        -(5*G*t*xe2^2+5*G*t*ye2^2-5*G*t*xe1*xe2+5*G*t*xe1*xe3-5*G*t*xe2*xe3-5*G*t*ye1*ye2+5*G*t*ye1*ye3-5*G*t*ye2*ye3)/(12*D)  (5*G*t*(xe1-xe2))/36  (5*G*t*(ye1-ye2))/36  (5*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(12*D)  (5*G*t*(xe1-xe2))/36  (5*G*t*(ye1-ye2))/36  (5*G*t*(xe1^2-2*xe1*xe2+xe2^2+ye1^2-2*ye1*ye2+ye2^2))/(12*D)  (5*G*t*(xe1-xe2))/36  (5*G*t*(ye1-ye2))/36  -(5*G*t*(xe1^2-2*xe1*xe2+xe2^2+ye1^2-2*ye1*ye2+ye2^2))/(9*D)  -(5*G*t*(xe1*xe2-xe1*xe3+xe2*xe3+ye1*ye2-ye1*ye3+ye2*ye3-xe2^2-ye2^2))/(9*D)  -(5*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(9*D)
        (5*G*t*(xe2-xe3))/36  (5*G*t*D)/144  0.0  -(5*G*t*(xe1-xe3))/36  (5*G*t*D)/144  0.0  (5*G*t*(xe1-xe2))/36  (5*G*t*D)/72  0.0  -(5*G*t*(xe1-xe2))/36  -(5*G*t*(xe1+xe2-2*xe3))/36  (5*G*t*(xe1+xe2-2*xe3))/36
        (5*G*t*(ye2-ye3))/36  0.0  (5*G*t*D)/144  -(5*G*t*(ye1-ye3))/36  0.0  (5*G*t*D)/144  (5*G*t*(ye1-ye2))/36  0.0  (5*G*t*D)/72  -(5*G*t*(ye1-ye2))/36  -(5*G*t*(ye1+ye2-2*ye3))/36  (5*G*t*(ye1+ye2-2*ye3))/36
        -(5*G*t*(xe1*xe2-xe1*xe3+xe2*xe3+ye1*ye2-ye1*ye3+ye2*ye3-xe2^2-ye2^2))/(9*D)  (5*G*t*(xe2-2*xe1+xe3))/36  (5*G*t*(ye2-2*ye1+ye3))/36  -(5*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(9*D)  -(5*G*t*(xe1-2*xe2+xe3))/36  -(5*G*t*(ye1-2*ye2+ye3))/36  -(5*G*t*(xe1^2-2*xe1*xe2+xe2^2+ye1^2-2*ye1*ye2+ye2^2))/(9*D)  -(5*G*t*(xe1-xe2))/36  -(5*G*t*(ye1-ye2))/36  -(10*G*t*(-xe1^2+xe1*xe2+xe1*xe3-xe2^2+xe2*xe3-xe3^2-ye1^2+ye1*ye2+ye1*ye3-ye2^2+ye2*ye3-ye3^2))/(9*D)  (10*G*t*(xe1*xe2-xe1*xe3+xe2*xe3+ye1*ye2-ye1*ye3+ye2*ye3-xe2^2-ye2^2))/(9*D)  (10*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(9*D)
        -(5*G*t*(xe2^2-2*xe2*xe3+xe3^2+ye2^2-2*ye2*ye3+ye3^2))/(9*D)  -(5*G*t*(xe2-xe3))/36  -(5*G*t*(ye2-ye3))/36  (5*G*t*(xe1*xe2-xe1*xe3-xe2*xe3+ye1*ye2-ye1*ye3-ye2*ye3+xe3^2+ye3^2))/(9*D)  (5*G*t*(xe1-2*xe2+xe3))/36  (5*G*t*(ye1-2*ye2+ye3))/36  -(5*G*t*(xe1*xe2-xe1*xe3+xe2*xe3+ye1*ye2-ye1*ye3+ye2*ye3-xe2^2-ye2^2))/(9*D)  -(5*G*t*(xe1+xe2-2*xe3))/36  -(5*G*t*(ye1+ye2-2*ye3))/36  (10*G*t*(xe1*xe2-xe1*xe3+xe2*xe3+ye1*ye2-ye1*ye3+ye2*ye3-xe2^2-ye2^2))/(9*D)  -(10*G*t*(-xe1^2+xe1*xe2+xe1*xe3-xe2^2+xe2*xe3-xe3^2-ye1^2+ye1*ye2+ye1*ye3-ye2^2+ye2*ye3-ye3^2))/(9*D)  -(10*G*t*(xe1*xe2-xe1*xe3-xe2*xe3+ye1*ye2-ye1*ye3-ye2*ye3+xe3^2+ye3^2))/(9*D)
        (5*G*t*(xe1*xe2-xe1*xe3-xe2*xe3+ye1*ye2-ye1*ye3-ye2*ye3+xe3^2+ye3^2))/(9*D)  -(5*G*t*(xe2-2*xe1+xe3))/36  -(5*G*t*(ye2-2*ye1+ye3))/36  -(5*G*t*(xe1^2-2*xe1*xe3+xe3^2+ye1^2-2*ye1*ye3+ye3^2))/(9*D)  (5*G*t*(xe1-xe3))/36  (5*G*t*(ye1-ye3))/36  -(5*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(9*D)  (5*G*t*(xe1+xe2-2*xe3))/36  (5*G*t*(ye1+ye2-2*ye3))/36  (10*G*t*(xe1*xe2+xe1*xe3-xe2*xe3+ye1*ye2+ye1*ye3-ye2*ye3-xe1^2-ye1^2))/(9*D)  -(10*G*t*(xe1*xe2-xe1*xe3-xe2*xe3+ye1*ye2-ye1*ye3-ye2*ye3+xe3^2+ye3^2))/(9*D)  -(10*G*t*(-xe1^2+xe1*xe2+xe1*xe3-xe2^2+xe2*xe3-xe3^2-ye1^2+ye1*ye2+ye1*ye3-ye2^2+ye2*ye3-ye3^2))/(9*D)
    ]

    # bending part (12×12) - only non-zero entries at rotation DOFs (rows/cols 2,3,5,6,8,9)
    kwtb = zeros(12, 12)
    # row 2 entries
    kwtb[2,2]=(t^3*(E22*xe2^2-2*E22*xe2*xe3+E22*xe3^2+G*ye2^2-2*G*ye2*ye3+G*ye3^2))/(24*D)
    kwtb[2,3]=(t^3*(xe2-xe3)*(ye2-ye3)*(E21+G))/(24*D)
    kwtb[2,5]=-(t^3*(E22*xe3^2+G*ye3^2+E22*xe1*xe2-E22*xe1*xe3-E22*xe2*xe3+G*ye1*ye2-G*ye1*ye3-G*ye2*ye3))/(24*D)
    kwtb[2,6]=-(t^3*(E21*xe2*ye1-E21*xe3*ye1-E21*xe2*ye3+E21*xe3*ye3+G*xe1*ye2-G*xe1*ye3-G*xe3*ye2+G*xe3*ye3))/(24*D)
    kwtb[2,8]=-(t^3*(E22*xe2^2+G*ye2^2-E22*xe1*xe2+E22*xe1*xe3-E22*xe2*xe3-G*ye1*ye2+G*ye1*ye3-G*ye2*ye3))/(24*D)
    kwtb[2,9]= (t^3*(E21*xe2*ye1-E21*xe2*ye2-E21*xe3*ye1+E21*xe3*ye2+G*xe1*ye2-G*xe1*ye3-G*xe2*ye2+G*xe2*ye3))/(24*D)
    # row 3
    kwtb[3,2]=kwtb[2,3]
    kwtb[3,3]=(t^3*(G*xe2^2-2*G*xe2*xe3+G*xe3^2+E11*ye2^2-2*E11*ye2*ye3+E11*ye3^2))/(24*D)
    kwtb[3,5]=-(t^3*(E12*xe1*ye2-E12*xe1*ye3-E12*xe3*ye2+E12*xe3*ye3+G*xe2*ye1-G*xe3*ye1-G*xe2*ye3+G*xe3*ye3))/(24*D)
    kwtb[3,6]=-(t^3*(E11*ye3^2+G*xe3^2+G*xe1*xe2-G*xe1*xe3-G*xe2*xe3+E11*ye1*ye2-E11*ye1*ye3-E11*ye2*ye3))/(24*D)
    kwtb[3,8]= (t^3*(E12*xe1*ye2-E12*xe1*ye3-E12*xe2*ye2+E12*xe2*ye3+G*xe2*ye1-G*xe2*ye2-G*xe3*ye1+G*xe3*ye2))/(24*D)
    kwtb[3,9]=-(t^3*(E11*ye2^2+G*xe2^2-G*xe1*xe2+G*xe1*xe3-G*xe2*xe3-E11*ye1*ye2+E11*ye1*ye3-E11*ye2*ye3))/(24*D)
    # row 5
    kwtb[5,2]=kwtb[2,5]; kwtb[5,3]=kwtb[3,5]
    kwtb[5,5]=(t^3*(E22*xe1^2-2*E22*xe1*xe3+E22*xe3^2+G*ye1^2-2*G*ye1*ye3+G*ye3^2))/(24*D)
    kwtb[5,6]=(t^3*(xe1-xe3)*(ye1-ye3)*(E21+G))/(24*D)
    kwtb[5,8]=-(t^3*(E22*xe1^2+G*ye1^2-E22*xe1*xe2-E22*xe1*xe3+E22*xe2*xe3-G*ye1*ye2-G*ye1*ye3+G*ye2*ye3))/(24*D)
    kwtb[5,9]=-(t^3*(E21*xe1*ye1-E21*xe1*ye2-E21*xe3*ye1+E21*xe3*ye2+G*xe1*ye1-G*xe2*ye1-G*xe1*ye3+G*xe2*ye3))/(24*D)
    # row 6
    kwtb[6,2]=kwtb[2,6]; kwtb[6,3]=kwtb[3,6]; kwtb[6,5]=kwtb[5,6]
    kwtb[6,6]=(t^3*(G*xe1^2-2*G*xe1*xe3+G*xe3^2+E11*ye1^2-2*E11*ye1*ye3+E11*ye3^2))/(24*D)
    kwtb[6,8]=-(t^3*(E12*xe1*ye1-E12*xe2*ye1-E12*xe1*ye3+E12*xe2*ye3+G*xe1*ye1-G*xe1*ye2-G*xe3*ye1+G*xe3*ye2))/(24*D)
    kwtb[6,9]=-(t^3*(E11*ye1^2+G*xe1^2-G*xe1*xe2-G*xe1*xe3+G*xe2*xe3-E11*ye1*ye2-E11*ye1*ye3+E11*ye2*ye3))/(24*D)
    # row 8
    kwtb[8,2]=kwtb[2,8]; kwtb[8,3]=kwtb[3,8]; kwtb[8,5]=kwtb[5,8]; kwtb[8,6]=kwtb[6,8]
    kwtb[8,8]=(t^3*(E22*xe1^2-2*E22*xe1*xe2+E22*xe2^2+G*ye1^2-2*G*ye1*ye2+G*ye2^2))/(24*D)
    kwtb[8,9]=(t^3*(xe1-xe2)*(ye1-ye2)*(E21+G))/(24*D)
    # row 9
    kwtb[9,2]=kwtb[2,9]; kwtb[9,3]=kwtb[3,9]; kwtb[9,5]=kwtb[5,9]; kwtb[9,6]=kwtb[6,9]; kwtb[9,8]=kwtb[8,9]
    kwtb[9,9]=(t^3*(G*xe1^2-2*G*xe1*xe2+G*xe2^2+E11*ye1^2-2*E11*ye1*ye2+E11*ye2^2))/(24*D)

    # static condensation
    inda = 1:9; indi = 10:12
    kwts_c = kwts[inda, inda] - kwts[inda, indi] * (kwts[indi, indi] \ kwts[indi, inda])
    kwtb_c = kwtb[inda, inda]

    # shear correction
    alpha = sum(diag(kwts_c[4:9, 4:9])) / max(sum(diag(kwtb_c[4:9, 4:9])), 1e-30)
    Cs = 0.0
    ke_wt = (1 / (1 + Cs*alpha)) * kwts_c + kwtb_c
    return ke_wt
end

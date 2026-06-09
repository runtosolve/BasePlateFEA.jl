function ke_uv_4n_condens_from_12to8dof_num(P1, P2, P3, P4, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp_in)
    xe1=P1[1]; ye1=P1[2]
    xe2=P2[1]; ye2=P2[2]
    xe3=P3[1]; ye3=P3[2]
    xe4=P4[1]; ye4=P4[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)
    D = [E11 E12 0.0; E21 E22 0.0; 0.0 0.0 G]
    Duv = D * t

    nGp = round(Int, sqrt(nGp_in))
    locx, locy, wei = GL4(nGp)

    kuv = zeros(12, 12)
    for i in 1:nGp, j in 1:nGp
        x = locx[i, j]; y = locy[i, j]
        J = [xe1*(y/4-1/4)-xe2*(y/4-1/4)+xe3*(y/4+1/4)-xe4*(y/4+1/4)  ye1*(y/4-1/4)-ye2*(y/4-1/4)+ye3*(y/4+1/4)-ye4*(y/4+1/4);
             xe1*(x/4-1/4)-xe2*(x/4+1/4)+xe3*(x/4+1/4)-xe4*(x/4-1/4)  ye1*(x/4-1/4)-ye2*(x/4+1/4)+ye3*(x/4+1/4)-ye4*(x/4-1/4)]
        detJ = det(J); Ji = inv(J)
        Ji11=Ji[1,1]; Ji12=Ji[1,2]; Ji21=Ji[2,1]; Ji22=Ji[2,2]
        N1x=(y-1)/4; N2x=-N1x; N3x=(y+1)/4; N4x=-N3x
        N5x=-2*x;   N6x=0.0
        N1y=(x-1)/4; N2y=-(x+1)/4; N3y=-N2y; N4y=-N1y
        N5y=0.0;    N6y=-2*y
        Buv = [
            N1x*Ji11+N1y*Ji12  0.0  N2x*Ji11+N2y*Ji12  0.0  N3x*Ji11+N3y*Ji12  0.0  N4x*Ji11+N4y*Ji12  0.0  N5x*Ji11+N5y*Ji12  0.0  N6x*Ji11+N6y*Ji12  0.0;
            0.0  N1x*Ji21+N1y*Ji22  0.0  N2x*Ji21+N2y*Ji22  0.0  N3x*Ji21+N3y*Ji22  0.0  N4x*Ji21+N4y*Ji22  0.0  N5x*Ji21+N5y*Ji22  0.0  N6x*Ji21+N6y*Ji22;
            N1x*Ji21+N1y*Ji22  N1x*Ji11+N1y*Ji12  N2x*Ji21+N2y*Ji22  N2x*Ji11+N2y*Ji12  N3x*Ji21+N3y*Ji22  N3x*Ji11+N3y*Ji12  N4x*Ji21+N4y*Ji22  N4x*Ji11+N4y*Ji12  N5x*Ji21+N5y*Ji22  N5x*Ji11+N5y*Ji12  N6x*Ji21+N6y*Ji22  N6x*Ji11+N6y*Ji12
        ]
        kuv += Buv' * Duv * Buv * detJ * wei[i, j]
    end

    inda = 1:8; indi = 9:12
    ke_uv = kuv[inda, inda] - kuv[inda, indi] * (kuv[indi, indi] \ kuv[indi, inda])
    return ke_uv
end

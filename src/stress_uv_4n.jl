function stress_uv_4n(d_elem, P1, P2, P3, P4, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp_in, ifgauss)
    xe1=P1[1]; ye1=P1[2]
    xe2=P2[1]; ye2=P2[2]
    xe3=P3[1]; ye3=P3[2]
    xe4=P4[1]; ye4=P4[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)

    if ifgauss
        nGp = round(Int, sqrt(nGp_in))
        locx, locy, wei = GL4(nGp)
    else
        nGp = 2
        locx = [-1.0  1.0;  1.0 -1.0]
        locy = [-1.0 -1.0;  1.0  1.0]
    end

    ii = 0
    str_uv = zeros(3, nGp * nGp)
    for i in 1:nGp, j in 1:nGp
        ii += 1
        x = locx[i, j]; y = locy[i, j]
        Jdet_expr = xe1*ye2-xe2*ye1-xe1*ye4+xe2*ye3-xe3*ye2+xe4*ye1+xe3*ye4-xe4*ye3
        Jdet_x    = -xe1*ye3+xe3*ye1+xe1*ye4+xe2*ye3-xe3*ye2-xe4*ye1-xe2*ye4+xe4*ye2
        Jdet_y    = -xe1*ye2+xe2*ye1+xe1*ye3-xe3*ye1-xe2*ye4+xe4*ye2+xe3*ye4-xe4*ye3
        Jdet_xy   = -xe1*ye3+xe3*ye1+xe1*ye4+xe2*ye3-xe3*ye2-xe4*ye1-xe2*ye4+xe4*ye2

        denom = Jdet_expr + x*Jdet_x + y*Jdet_y + x*y*(xe1*ye3-xe3*ye1-xe1*ye4-xe2*ye3+xe3*ye2+xe4*ye1+xe2*ye4-xe4*ye2+xe1*ye3-xe3*ye1-xe1*ye4-xe2*ye3+xe3*ye2+xe4*ye1+xe2*ye4-xe4*ye2)
        # Use the standard Jacobian formulation for 4-node element
        J11 = xe1*(y/4-1/4) - xe2*(y/4-1/4) + xe3*(y/4+1/4) - xe4*(y/4+1/4)
        J12 = ye1*(y/4-1/4) - ye2*(y/4-1/4) + ye3*(y/4+1/4) - ye4*(y/4+1/4)
        J21 = xe1*(x/4-1/4) - xe2*(x/4+1/4) + xe3*(x/4+1/4) - xe4*(x/4-1/4)
        J22 = ye1*(x/4-1/4) - ye2*(x/4+1/4) + ye3*(x/4+1/4) - ye4*(x/4-1/4)
        Jmat = [J11 J12; J21 J22]
        detJ = det(Jmat)
        Ji = inv(Jmat)
        Ji11=Ji[1,1]; Ji12=Ji[1,2]; Ji21=Ji[2,1]; Ji22=Ji[2,2]

        N1x=(y-1)/4; N2x=-(y-1)/4; N3x=(y+1)/4; N4x=-(y+1)/4
        N1y=(x-1)/4; N2y=-(x+1)/4; N3y=(x+1)/4; N4y=-(x-1)/4

        Buv = [
            N1x*Ji11+N1y*Ji12  0.0  N2x*Ji11+N2y*Ji12  0.0  N3x*Ji11+N3y*Ji12  0.0  N4x*Ji11+N4y*Ji12  0.0;
            0.0  N1x*Ji21+N1y*Ji22  0.0  N2x*Ji21+N2y*Ji22  0.0  N3x*Ji21+N3y*Ji22  0.0  N4x*Ji21+N4y*Ji22;
            N1x*Ji21+N1y*Ji22  N1x*Ji11+N1y*Ji12  N2x*Ji21+N2y*Ji22  N2x*Ji11+N2y*Ji12  N3x*Ji21+N3y*Ji22  N3x*Ji11+N3y*Ji12  N4x*Ji21+N4y*Ji22  N4x*Ji11+N4y*Ji12
        ]
        D_mat = [E11 E12 0.0; E21 E22 0.0; 0.0 0.0 G]
        str_uv[:, ii] = D_mat * Buv * d_elem / t
    end
    return str_uv
end

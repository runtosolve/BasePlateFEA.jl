function ke_wt_4n_condens_from_18to12dof_num(P1, P2, P3, P4, t, Ex, Ey, nuxy, nuyx, G, mat_sym, nGp_in)
    X1=P1[1]; Y1=P1[2]
    X2=P2[1]; Y2=P2[2]
    X3=P3[1]; Y3=P3[2]
    X4=P4[1]; Y4=P4[2]

    E11 = Ex  / (1 - nuxy*nuyx)
    E22 = Ey  / (1 - nuxy*nuyx)
    E12 = nuyx*Ex / (1 - nuxy*nuyx)
    E21 = mat_sym == 1 ? E12 : nuxy*Ey / (1 - nuxy*nuyx)
    D   = [E11 E12 0.0; E21 E22 0.0; 0.0 0.0 G]
    Dwtb = D * t^3 / 12
    Dwts = Matrix{Float64}(I, 2, 2) * (5/6 * G * t)

    nGp = round(Int, sqrt(nGp_in))
    locx, locy, wei = GL4(nGp)

    kwts = zeros(18, 18)
    kwtb = zeros(18, 18)

    for i in 1:nGp, j in 1:nGp
        x = locx[i, j]; y = locy[i, j]
        J = [X1*(y/4-1/4)-X2*(y/4-1/4)+X3*(y/4+1/4)-X4*(y/4+1/4)  Y1*(y/4-1/4)-Y2*(y/4-1/4)+Y3*(y/4+1/4)-Y4*(y/4+1/4);
             X1*(x/4-1/4)-X2*(x/4+1/4)+X3*(x/4+1/4)-X4*(x/4-1/4)  Y1*(x/4-1/4)-Y2*(x/4+1/4)+Y3*(x/4+1/4)-Y4*(x/4-1/4)]
        detJ = det(J); Ji = inv(J)
        Ji11=Ji[1,1]; Ji12=Ji[1,2]; Ji21=Ji[2,1]; Ji22=Ji[2,2]
        N1=(x-1)*(y-1)/4; N2=-(x+1)*(y-1)/4; N3=(x+1)*(y+1)/4; N4=-(x-1)*(y+1)/4
        N5=1-x^2; N6=1-y^2
        N1x=(y-1)/4; N2x=-N1x; N3x=(y+1)/4; N4x=-N3x; N5x=-2*x; N6x=0.0
        N1y=(x-1)/4; N2y=-(x+1)/4; N3y=-N2y; N4y=-N1y; N5y=0.0; N6y=-2*y

        Bwts = [
            N1x*Ji11+N1y*Ji12  0.0  N1  N2x*Ji11+N2y*Ji12  0.0  N2  N3x*Ji11+N3y*Ji12  0.0  N3  N4x*Ji11+N4y*Ji12  0.0  N4  N5x*Ji11+N5y*Ji12  0.0  N5  N6x*Ji11+N6y*Ji12  0.0  N6;
            N1x*Ji21+N1y*Ji22  -N1  0.0  N2x*Ji21+N2y*Ji22  -N2  0.0  N3x*Ji21+N3y*Ji22  -N3  0.0  N4x*Ji21+N4y*Ji22  -N4  0.0  N5x*Ji21+N5y*Ji22  -N5  0.0  N6x*Ji21+N6y*Ji22  -N6  0.0
        ]
        kwts += Bwts' * Dwts * Bwts * detJ * wei[i, j]

        Bwtb = zeros(3, 18)
        # row 1: d(ty)/dx
        for k in 1:6
            Nkx = [N1x, N2x, N3x, N4x, N5x, N6x][k]
            Nky = [N1y, N2y, N3y, N4y, N5y, N6y][k]
            col = (k-1)*3 + 2
            Bwtb[1, col] = Nkx*Ji11 + Nky*Ji12
        end
        # row 2: -d(tx)/dy  (note tx is 3rd DOF per node)
        for k in 1:6
            Nkx = [N1x, N2x, N3x, N4x, N5x, N6x][k]
            Nky = [N1y, N2y, N3y, N4y, N5y, N6y][k]
            col = (k-1)*3 + 3
            Bwtb[2, col] = -(Nkx*Ji21 + Nky*Ji22)
        end
        # row 3: d(ty)/dy - d(tx)/dx  -> Bwtb[3,:] = -Bwtb[1,tx_cols] and Bwtb[3,ty_cols] refers to col+1
        for k in 1:6
            Nkx = [N1x, N2x, N3x, N4x, N5x, N6x][k]
            Nky = [N1y, N2y, N3y, N4y, N5y, N6y][k]
            col_ty = (k-1)*3 + 2
            col_tx = (k-1)*3 + 3
            Bwtb[3, col_ty] = Nkx*Ji21 + Nky*Ji22
            Bwtb[3, col_tx] = -(Nkx*Ji11 + Nky*Ji12)
        end
        kwtb += Bwtb' * Dwtb * Bwtb * detJ * wei[i, j]
    end

    inda = 1:12; indi = 13:18
    kwt1 = kwts + kwtb
    ke_wt = kwt1[inda, inda] - kwt1[inda, indi] * (kwt1[indi, indi] \ kwt1[indi, inda])
    return ke_wt
end

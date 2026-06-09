function lv_el4_edge(elems_row, nodes, px, py, pz, dpn)
    T, P1e, P2e, P3e, P4e = ct_4node_g2e(elems_row, nodes)
    xe1=P1e[1]; ye1=P1e[2]
    xe2=P2e[1]; ye2=P2e[2]
    xe3=P3e[1]; ye3=P3e[2]
    xe4=P4e[1]; ye4=P4e[2]

    L12 = sqrt((xe1-xe2)^2 + (ye1-ye2)^2)
    L23 = sqrt((xe2-xe3)^2 + (ye2-ye3)^2)
    L34 = sqrt((xe3-xe4)^2 + (ye3-ye4)^2)
    L14 = sqrt((xe1-xe4)^2 + (ye1-ye4)^2)

    lv0_all = [
        [L12/3, L12/6, 0.0, 0.0],
        [L12/6, L12/3, 0.0, 0.0],
        [0.0, L23/3, L23/6, 0.0],
        [0.0, L23/6, L23/3, 0.0],
        [0.0, 0.0, L34/3, L34/6],
        [0.0, 0.0, L34/6, L34/3],
        [L14/3, 0.0, 0.0, L14/6],
        [L14/6, 0.0, 0.0, L14/3],
    ]

    elemnode = 4
    lv = zeros(elemnode * dpn)
    for i in 1:2*elemnode
        lv0 = lv0_all[i]
        pe = zeros(dpn)
        pe[1:min(3,dpn)] = [px[i], py[i], pz[i]][1:min(3,dpn)]
        pp = zeros(dpn*elemnode, elemnode)
        for k in 1:elemnode
            pp[(k-1)*dpn+1:k*dpn, k] = pe
        end
        lv += pp * lv0
    end
    return lv
end

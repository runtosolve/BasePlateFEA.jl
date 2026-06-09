function FEanal_incr(iincr, dincr, d_ini, Ke, lv, nodes, el_nodes, el_props,
                     pl_el_edges, pl_edge_nodes, nodenr,
                     PLATE_table, MAT_table, SUPPORT_table, CONDSUP_table, LOAD_table,
                     FE_dat, mat_sym, dpn, tolers)
    if iincr == 1
        println("Stiffness matrix calculation ")
        t0 = time()
        Ke = stiffmat_e(nodes, el_nodes, el_props, nodenr, MAT_table, FE_dat, mat_sym, dpn)
        println("  done in $(round(time()-t0, digits=2)) s")

        println("Adding supports ")
        t0 = time()
        Ke_sup = fem_support(Ke, PLATE_table, SUPPORT_table, FE_dat, nodes, el_nodes, el_props, nodenr, dpn, tolers)
        Ke = Ke + Ke_sup
        println("  done in $(round(time()-t0, digits=2)) s")

        println("Load vector ")
        t0 = time()
        lv = load_vec_UDL(nodes, el_nodes, el_props, nodenr, PLATE_table, LOAD_table, FE_dat, dpn)
        lv .+= load_vec_edge(nodes, el_nodes, el_props, pl_el_edges, pl_edge_nodes, nodenr, PLATE_table, LOAD_table, FE_dat, dpn, tolers)
        lv .+= load_vec_point(nodes, nodenr, PLATE_table, LOAD_table, dpn, tolers)
        println("  done in $(round(time()-t0, digits=2)) s")
    end

    println("Load increment $iincr")
    t0 = time()
    dlv = lv * dincr[iincr]

    if iincr == 1
        uncond = true
        Ke_cond = cond_support(Ke, PLATE_table, CONDSUP_table, FE_dat, nodes, el_nodes, el_props, nodenr, d_ini, dpn, uncond, tolers)
        d_ini = Array(Ke + Ke_cond) \ dlv
    end

    uncond = false
    Ke_cond = cond_support(Ke, PLATE_table, CONDSUP_table, FE_dat, nodes, el_nodes, el_props, nodenr, d_ini, dpn, uncond, tolers)
    dd = Array(Ke + Ke_cond) \ dlv
    if iincr == 1; d_ini = zeros(length(d_ini)); end
    d_stat = d_ini + dd
    println("  done in $(round(time()-t0, digits=2)) s")

    return d_stat, Ke, lv
end

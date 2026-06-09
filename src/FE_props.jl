function FE_props()
    fe_dat = Dict{Int, FEData}()

    # t3
    fe_dat[11] = FEData(3, 3, 6, 6, false, false, 0, "Mindlin", 12, 9, true, false, 0, 100, 18, 1)

    # t31
    fe_dat[12] = FEData(3, 3, 8, 6, true, false, 0, "Mindlin", 12, 9, true, false, 0, 100, 18, 1)

    # q4
    fe_dat[21] = FEData(4, 4, 8, 8, false, true, 4, "Mindlin", 18, 12, true, true, 9, 100, 24, 4)

    # q42
    fe_dat[22] = FEData(4, 4, 12, 8, true, true, 4, "Mindlin", 18, 12, true, true, 9, 100, 24, 4)

    return fe_dat
end

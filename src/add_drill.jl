function add_drill(k5dpn::AbstractMatrix, nodenr::Int, ifelastic::Bool)
    ind = Int[]
    i0 = 0
    for i in 1:nodenr
        append!(ind, (1:5) .+ i0)
        i0 += 6
    end
    k6dpn = zeros(6 * nodenr, 6 * nodenr)
    k6dpn[ind, ind] = k5dpn

    if ifelastic
        ind2 = Int[]
        i0 = 0
        for i in 1:nodenr
            append!(ind2, [4, 5] .+ i0)
            i0 += 5
        end
        kd = diag(k5dpn)
        stif = minimum(kd[ind2]) / 100
        for i in 1:nodenr
            k6dpn[6*i, 6*i] = stif
        end
    end
    return k6dpn
end

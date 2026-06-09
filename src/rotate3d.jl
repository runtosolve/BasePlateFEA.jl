function rotate3d(T3::AbstractMatrix, nodenr::Int, dpn::Int)
    T = Matrix{Float64}(I, nodenr * dpn, nodenr * dpn)
    if dpn == 5
        T2 = T3[1:2, 1:2]
        i0 = 0
        for i in 1:nodenr
            ind = (1:3) .+ i0;  T[ind, ind] = T3
            ind = (4:5) .+ i0;  T[ind, ind] = T2
            i0 += dpn
        end
    elseif dpn == 6
        i0 = 0
        for i in 1:nodenr
            ind = (1:3) .+ i0;  T[ind, ind] = T3
            ind = (4:6) .+ i0;  T[ind, ind] = T3
            i0 += dpn
        end
    end
    return T
end

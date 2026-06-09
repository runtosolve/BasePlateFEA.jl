function mesh4_rect(P1, P2, P3, P4, p_anchor_mat, size_x, size_y, tol)
    nmin = 1
    P0, T = plate_transform(P1, P2, P3, P4)

    P1s = T' * (P1 - P0)
    P2s = T' * (P2 - P0)
    P3s = T' * (P3 - P0)
    P4s = T' * (P4 - P0)

    x_anch = [P1s[1], P3s[1]]
    y_anch = [P1s[2], P3s[2]]
    for i in 1:size(p_anchor_mat, 1)
        Ps = T' * (p_anchor_mat[i, :] - P0)
        if minimum(abs.(x_anch .- Ps[1])) > tol
            push!(x_anch, Ps[1])
        end
        if minimum(abs.(y_anch .- Ps[2])) > tol
            push!(y_anch, Ps[2])
        end
    end
    x_anch = sort(x_anch)
    y_anch = sort(y_anch)

    x_seed = [x_anch[1]]
    for i in 2:length(x_anch)
        nx = max(round(Int, (x_anch[i] - x_anch[i-1]) / size_x), nmin)
        dx = (x_anch[i] - x_anch[i-1]) / nx
        append!(x_seed, x_seed[end] .+ (1:nx) .* dx)
    end
    nrx = length(x_seed)

    y_seed = [y_anch[1]]
    for i in 2:length(y_anch)
        ny = max(round(Int, (y_anch[i] - y_anch[i-1]) / size_y), nmin)
        dy = (y_anch[i] - y_anch[i-1]) / ny
        append!(y_seed, y_seed[end] .+ (1:ny) .* dy)
    end
    nry = length(y_seed)

    a = x_seed[nrx] - x_seed[1]
    b = y_seed[nry] - y_seed[1]
    x0 = x_seed[1]; y0 = y_seed[1]

    nodes = zeros(nrx * nry, 3)
    nr = 0
    for j in 1:nry, i in 1:nrx
        x = (x_seed[i] - x0) / a
        y = (y_seed[j] - y0) / b
        w1 = (1-x)*(1-y); w2 = x*(1-y); w3 = x*y; w4 = (1-x)*y
        Ps = w1*P1s[1:3] + w2*P2s[1:3] + w3*P3s[1:3] + w4*P4s[1:3]
        nr += 1
        nodes[nr, :] = T * Ps + P0
    end

    nelx = nrx - 1; nely = nry - 1
    elems = zeros(Int, nelx * nely, 4)
    nr = 0
    for j in 1:nely, i in 1:nelx
        nr += 1
        elems[nr, 1] = nrx*(j-1) + i
        elems[nr, 2] = nrx*(j-1) + i + 1
        elems[nr, 3] = nrx*j     + i + 1
        elems[nr, 4] = nrx*j     + i
    end
    return nodes, elems
end

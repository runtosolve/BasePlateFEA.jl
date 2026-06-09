function GL4(np::Int)
    if np == 1
        locx = reshape([0.0], 1, 1)
        locy = reshape([0.0], 1, 1)
        wei  = reshape([4.0], 1, 1)
    elseif np == 2
        g = 1/sqrt(3)
        locx = [-g  -g;  g   g]
        locy = [-g   g; -g   g]
        wei  = [ 1.0  1.0;  1.0  1.0]
    elseif np == 3
        g = sqrt(3/5)
        locx = [-g  -g  -g;  0   0   0;  g   g   g]
        locy = [-g   0   g; -g   0   g; -g   0   g]
        w1 = [5.0, 8.0, 5.0]
        w2 = [5.0, 8.0, 5.0]
        wei = w1 * w2' / 81
    elseif np == 4
        g1 = sqrt(3/7 + 2/7*sqrt(6/5))
        g2 = sqrt(3/7 - 2/7*sqrt(6/5))
        locx = [-g1  -g1  -g1  -g1;
                -g2  -g2  -g2  -g2;
                 g2   g2   g2   g2;
                 g1   g1   g1   g1]
        locy = [-g1  -g2   g2   g1;
                -g1  -g2   g2   g1;
                -g1  -g2   g2   g1;
                -g1  -g2   g2   g1]
        wv = [18-sqrt(30), 18+sqrt(30), 18+sqrt(30), 18-sqrt(30)]
        wei = (wv * wv') / 36^2
    else
        error("GL4: unsupported np=$np")
    end
    return locx, locy, wei
end

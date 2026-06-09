function GL3(np::Int)
    if np == 1
        locx = [1/3]
        locy = [1/3]
        wei  = [1/2]
    elseif np == 4
        locx = [0.1889958, 0.7053418, 0.1279915, 0.4776709]
        locy = [0.1889958, 0.1279915, 0.7053418, 0.4776709]
        wei  = [0.1971688, 0.125,     0.125,     0.05283122]
    else
        error("GL3: unsupported np=$np")
    end
    return locx, locy, wei
end

function cs_def(h, b, c, tb, th, tc)
    cs_node = [
        1.0  b      -h/2+c
        2.0  b      -h/2
        3.0  0.0    -h/2
        4.0  0.0     h/2
        5.0  b       h/2
        6.0  b       h/2-c
    ]
    cs_elem = [
        1  1  2  tc
        2  2  3  tb
        3  3  4  th
        4  4  5  tb
        5  5  6  tc
    ]
    return cs_node, Float64.(cs_elem)
end

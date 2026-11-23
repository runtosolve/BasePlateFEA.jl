# --------------------------------------------------------------------------------------------
# Thin-walled section properties of the lipped-C upright (cs_def.m, cs_prop.m, cs_data.m).
# --------------------------------------------------------------------------------------------

"""
    cs_data(h, b, c, t) -> (A, yCG)

Area and distance of the centroid from the web of a lipped C with centerline web depth `h`,
flange width `b`, lip `c` and thickness `t` (straight-line model), as `cs_data.m`.
"""
function cs_data(h, b, c, t)
    A = t * (h + 2b + 2c)
    return A, t * (2b * b / 2 + 2c * b) / A
end

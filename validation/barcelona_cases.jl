# Barcelona base plate tests (Roure et al.): model definitions used for the validation runs.
# Geometry follows S. Ádány's bp_Barcelona.m; thicknesses and steel strengths are the measured
# values in "Bases 2011-2012_material properties and other test data from Francesc Roure.pdf".
# The strain at the tensile strength is not reported: 15 % (plate) and 10 % (upright) are assumed.
using BasePlateFEA

const BARCELONA = Dict(
    :A => (t = 0.363, tu = 0.0984, plate = (61.839, 90.547, 0.15), upright = (64.450, 71.802, 0.10),
           loads = [5.0, 10.0, 15.0], θmax = 0.08),
    :B => (t = 0.372, tu = 0.1016, plate = (47.891, 73.424, 0.15), upright = (68.911, 77.810, 0.10),
           loads = [5.0, 7.5, 10.0, 12.5, 15.0, 20.0], θmax = 0.08),
    :C => (t = 0.1291, tu = 0.1016, plate = (53.226, 64.273, 0.15), upright = (68.911, 77.810, 0.10),
           loads = [5.0, 7.5, 10.0, 12.5, 15.0, 20.0], θmax = 0.05),
)

"geometry of `bp_Barcelona.m` with the measured thicknesses"
function barcelona_test_geometry(type::Symbol)
    g = barcelona_geometry(type); c = BARCELONA[type]
    return BasePlateGeometry(W = g.W, L = g.L, t = c.t, Da = g.Da, Dc = g.Dc, Sa = g.Sa, R = g.R,
                             h = g.h, b = g.b, c = g.c, tu = c.tu)
end

"model input: `plastic = false` is the elastic model with the settings of `bp_Barcelona.m`"
function barcelona_input(type::Symbol; plastic::Bool, kwargs...)
    c = BARCELONA[type]
    plastic || return BasePlateInput(; geometry = barcelona_test_geometry(type), kwargs...)
    return BasePlateInput(; geometry = barcelona_test_geometry(type), plate_steel = bilinear_steel(c.plate...),
                          upright_steel = bilinear_steel(c.upright...), kwargs...)
end

loadtag(N) = replace(string(N), r"\.0$" => "", "." => "p")

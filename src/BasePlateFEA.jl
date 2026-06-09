module BasePlateFEA

using LinearAlgebra
using SparseArrays
using Statistics
using Random
using GLMakie
using GeometryBasics: Point3f, TriangleFace, Mesh

include("types.jl")
include("FE_props.jl")
include("GL3.jl")
include("GL4.jl")
include("plate_transform.jl")
include("ct_3node_g2e.jl")
include("ct_4node_g2e.jl")
include("rotate3d.jl")
include("add_drill.jl")
include("cs_def.jl")
include("cs_prop.jl")
include("cs_data.jl")
include("ke_uv_3n_nocondens_6dof.jl")
include("ke_uv_3n_condens_from_8to6dof.jl")
include("ke_uv_4n_nocondens_8dof_num.jl")
include("ke_uv_4n_condens_from_12to8dof_num.jl")
include("ke_wt_3n_condens_from_12to9dof_MIN.jl")
include("ke_wt_4n_condens_from_18to12dof_num.jl")
include("stress_uv_3n.jl")
include("stress_uv_4n.jl")
include("lv_el3_edge.jl")
include("lv_el4_edge.jl")
include("collect_anchor.jl")
include("mesh4_rect.jl")
include("edges_by_plate.jl")
include("elem_divide_by_plate.jl")
include("four2three_by_plate.jl")
include("stiffmat_e.jl")
include("fem_support.jl")
include("cond_support.jl")
include("load_vec_UDL.jl")
include("load_vec_edge.jl")
include("load_vec_point.jl")
include("main_modelgen.jl")
include("main_stressout.jl")
include("FEanal_incr.jl")
include("plot_res_tri.jl")

export Tolerances, FEData, HoleData
export FE_props
export GL3, GL4
export plate_transform, ct_3node_g2e, ct_4node_g2e, rotate3d, add_drill
export cs_def, cs_prop, cs_data
export ke_uv_3n_nocondens_6dof, ke_uv_3n_condens_from_8to6dof
export ke_uv_4n_nocondens_8dof_num, ke_uv_4n_condens_from_12to8dof_num
export ke_wt_3n_condens_from_12to9dof_MIN, ke_wt_4n_condens_from_18to12dof_num
export stress_uv_3n, stress_uv_4n
export lv_el3_edge, lv_el4_edge
export collect_anchor, mesh4_rect, edges_by_plate
export elem_divide_by_plate, four2three_by_plate, HOLE_TABLE_entry
export stiffmat_e, fem_support, cond_support
export load_vec_UDL, load_vec_edge, load_vec_point
export main_modelgen, main_stressout, FEanal_incr
export plot_res_tri

end  # module 

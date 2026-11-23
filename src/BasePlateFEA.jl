"""
    BasePlateFEA

Shell finite element model of a cold-formed steel rack upright on a base plate, anchored to a
concrete block, for the moment–rotation response of the column base. Julia implementation of
S. Ádány's MATLAB column-base model (`bp_Barcelona.m`, November 2026): DistMesh triangles on the
base plate and quadrilaterals on the upright, rigid links at the loaded cross-section and at the
bolt holes, an elastic half-space under compression-only contact links with a yield stress,
tension-only anchor bolt springs, and a two-step (axial force, then moment) incremental analysis.

The shell elements are those of TriShellFiniteElement.jl and QuadShellFiniteElement.jl; their
geometrically nonlinear, layered von Mises plasticity formulation is used for the parts that are
given a plastic steel.
"""
module BasePlateFEA

using LinearAlgebra, SparseArrays, Printf
import Random
import StressStrainCurves
using DelimitedFiles: readdlm
using StaticArrays
using Ferrite
using Tensors: Vec
import DelaunayTriangulation
import JSON3
import TriShellFiniteElement as TSF
import QuadShellFiniteElement as QSF

export BasePlateGeometry, BasePlateInput, Steel, bilinear_steel, yun_gardner_steel, SpringLaw
export barcelona_geometry, barcelona_Mmax0
export Mesh, generate_mesh, mesh_from_matlab, apply_imperfection, buckling_modes
export Model, analyze, Result, Frame, moment_rotation, base_moment
export element_fields, nodal_average, viewer_data, write_viewer

include("input.jl")
include("section.jl")
include("distmesh.jl")
include("mesh.jl")
include("parts.jl")
include("shell.jl")
include("solve.jl")
include("buckling.jl")
include("results.jl")
include("viewer.jl")

end # module BasePlateFEA

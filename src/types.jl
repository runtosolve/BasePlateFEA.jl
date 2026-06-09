"""
Types used throughout the OneRack FE framework.
"""

mutable struct Tolerances
    zero_thickn::Float64
    zero_spring::Float64
    zero_dist::Float64
    zeroFElength::Float64
end

mutable struct FEData
    shape::Int
    nodenr::Int
    m_dof_ini::Int
    m_dof_final::Int
    m_ifcondens::Bool
    m_numint::Bool
    m_nGpe::Int
    b_theory::String
    b_dof_ini::Int
    b_dof_final::Int
    b_ifcondens::Bool
    b_numint::Bool
    b_nGpe::Int
    nonzer_ke::Int
    nonzer_lv::Int
    nGpg::Int
end

mutable struct HoleData
    coord::Vector{Vector{Float64}}
    radius::Vector{Float64}
    sector::Vector{Int}
    holes::Vector{Matrix{Float64}}
end

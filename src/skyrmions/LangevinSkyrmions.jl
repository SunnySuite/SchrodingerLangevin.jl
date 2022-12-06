module LangevinSkyrmions
 
using Sunny
using LinearAlgebra, Random


################################################################################
# Function to create model in Sunny 
################################################################################
"""
    su3_skyrmion_model(dims;
    J_ref=1.0, h = 15.5, D = 19.0, lat_type = 1, rng = nothing,
)

Creates a Sunny implementation of the skyrmion-supporting spin-1 model described in
Zhang et al. (2022). 
"""
function su3_skyrmion_model(dims;
    J_ref=1.0, h = 15.5, D = 19.0, lat_type = 1, rng = nothing,
)
    isnothing(rng) && (rng = MersenneTwister())
    h *= 8.637992737026183 # × (1/2*μᵦ) to make unitless

    #= Model parameters =#
    J₁ = -J_ref
    J₂ = (2.0/(1+√5)) * J_ref 
    Δ = 2.6
    N = 3

    #= Build crystal =#
    α, β, γ = 90, 90, 60 
    a, b, c = 1.0, 1.0, 2.0
    lat_vecs = Sunny.lattice_vectors(a, b, c, α, β, γ)
    basis_vecs = [[0,0,0]]
    cryst = Crystal(lat_vecs, basis_vecs)


    #= Set up Hamiltonian =#
    # Exchange
    bond1 = Bond(1, 1, [1, 0, 0])  #  nearest
    bond2 = Bond(1, 1, [-1, 2, 0])  # next-nearest
    Jex1 = J₁ * [1.0 0.0 0.0;
                 0.0 1.0 0.0;
                 0.0 0.0 Δ]
    Jex2 = J₂ * [1.0 0.0 0.0;
                 0.0 1.0 0.0;
                 0.0 0.0 Δ]
    ex1 = exchange(Jex1, bond1)
    ex2 = exchange(Jex2, bond2)

    # External field
    field = external_field([0.0 0.0 h])

    # Anisotropy
    Λ = D * 𝒮[3]^2
    aniso = anisotropy(Λ, 1)

    # Spin system
    interactions = [ex1, ex2, field, aniso]
    site_infos = [SiteInfo(1; N=3)]
    
    return SpinSystem(cryst, interactions, dims, site_infos; rng)
end


################################################################################
# Functions to generate figures and animations 
################################################################################

using GLMakie, Observables, ColorSchemes, ColorTypes

export fig3, quench_animation, thermalize_animation

include("plot_utils.jl")
include("fig3.jl")
include("animations.jl")


end
using Random

function su3_skyrmion_model(dims;
    J_ref=1.0,
    h = 15.5,
    D = 19.0,
    rescaling = 1.0,
    rng = nothing,
)
    isnothing(rng) && (rng = MersenneTwister())
    h *= 8.637992737026183 # Make field unitless × (1/2*μᵦ)

    #= Model parameters =#
    J₁ = -J_ref
    J₂ = (2.0/(1+√(5*rescaling))) * J_ref 
    # J₂ = 0.618034 
    Δ = 2.6
    N = 3


    #= Build crystal =#
    # α, β, γ = 90, 90, 60 
    α, β, γ = 90, 90, 120
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
    S = Sunny.gen_spin_ops(N)
    Λ = D * (S[3]^2)
    aniso = SUN_anisotropy(Λ, 1)

    # Spin system
    interactions = [ex1, ex2, field, aniso]
    site_infos = [SiteInfo(1; N=3)]
    
    return SpinSystem(cryst, interactions, dims, site_infos; rng)
end
using Test
using SchrodingerMidpoint
using LinearAlgebra


# Levi-Civita symbol
const ϵ = [(i-j)*(j-k)*(k-i)/2 for i=1:3, j=1:3, k=1:3]

# Kronecker delta
δ(i,j) = (i==j) ? 1 : 0


### Verify 𝔰𝔲(2) irreps

@testset "su(2) properties" begin

    for N = 2:5
        S₀ = (N-1)/2
        S = spin_operators(N)

        for i=1:3, j=1:3
            # Test commutation relations
            @test S[i]*S[j] - S[j]*S[i] ≈ im * sum(ϵ[i,j,k]*S[k] for k=1:3)

            # Test orthonormality
            @test tr(S[i]*S[j]) ≈ (2/3)*S₀*(S₀+1/2)*(S₀+1)*δ(i,j)
        end

        # Test magnitude
        @test sum(S[i]^2 for i=1:3) ≈ S₀*(S₀+1)*I

        # Test dipole -> ket -> dipole round trip
        n = S₀ * normalize(randn(Vec3))
        ψ = SchrodingerMidpoint.infer_ket_from_dipole(S, n)
        @test SchrodingerMidpoint.spin_bilinear(S, ψ) ≈ n
    end


    # Interesting property of spin-1/2 and spin-1 operators

    a = randn(Vec3{Float64})
    b = randn(Vec3{Float64})

    N = 2
    S = spin_operators(N)
    @test (a⋅S)*(b⋅S)*(a⋅S) ≈ (a⋅b)*(a⋅S)/2 - (a⋅a)*(b⋅S)/4 

    N = 3
    S = spin_operators(N)
    @test (a⋅S)*(b⋅S)*(a⋅S) ≈ (a⋅b)*(a⋅S)

end


### Approximate energy conservation

@testset "Energy conservation" begin
    function build_sys(; N, L)
        # Arbitrary site-dependent anisotropies
        S = spin_operators(N)
        Λ = [cos(θ) * (S[1]*S[2]^2 + S[2]^2*S[1]) for θ = range(0, 2π, length=L)]
    
        # Some complicated spin chain
        sys = System(; N, L, J=1., B=[0, 0, 1], D=1., Λ)
    
        # Modify spin at i=1
        θ = π/4
        set_state_from_dipole!(sys, 1, Vec3(sin(θ), 0, cos(θ)))

        return sys
    end

    # Symplectic integrator avoids energy drift over arbitrarily long
    # time-scales
    N = 4
    L = 5
    Δt = 0.1
    sys = build_sys(; N, L)
    E₀ = energy(sys)
    for i = 1:1000
        midpoint_step!(sys, Δt)
    end
    @test isapprox(E₀,  energy(sys); rtol=0.05)

    # Non-symplectic integrator needs a precisely controlled trajectory length
    # for this test
    N = 4
    L = 5
    Δt = 0.02
    sys = build_sys(; N, L)
    E₀ = energy(sys)
    for i = 1:1000
        heun_step!(sys, Δt)
    end
    @test isapprox(E₀,  energy(sys); rtol=0.02)

    # LLD case: dipole only
    N = 0
    L = 5
    Δt = 0.1
    sys = build_sys(; N, L)
    E₀ = energy(sys)
    for i = 1:1000
        spherical_midpoint_step!(sys, Δt)
    end
    @test isapprox(E₀,  energy(sys); rtol=0.01)
end

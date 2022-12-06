module SchrodingerLangevin

export Vec3, Mat3, spin_operators, System, set_state_from_dipole!, energy, heun_langevin_step!, rand!

using LinearAlgebra
using StaticArrays
using Random


const Vec3 = SVector{3, T} where T
const Mat3 = SMatrix{3, 3, T} where T

# Dot product on 3-vectors of operators
function LinearAlgebra.dot(a::Vec3{Float64}, b::Vec3{T}) where T
    return a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
end

# Generators for 𝔰𝔲(2) in the irreducible, N-dimensional representation
function spin_operators(N)
    if iszero(N)
        Sx = Sy = Sz = SMatrix{0, 0, ComplexF64, 0}()
    else
        S₀ = (N-1)/2
        a = 1:N-1
        off = @. sqrt(2(S₀+1)*a - a*(a+1)) / 2
        Sx = SMatrix{N, N, ComplexF64, N^2}(diagm(1 => off, -1 => off))
        Sy = SMatrix{N, N, ComplexF64, N^2}(diagm(1 => -im*off, -1 => +im*off))
        Sz = SMatrix{N, N, ComplexF64, N^2}(diagm((N-1)/2 .- (0:N-1)))
    end
    return Vec3(Sx, Sy, Sz)
end

function spin_bilinear(S, Z)
    return real.(Vec3(Z'*S[1]*Z, Z'*S[2]*Z, Z'*S[3]*Z))
end

function infer_ket_from_dipole(S, n)
    # Find a ket (up to an irrelevant phase) that corresponds to a pure dipole.
    # TODO, we can do this much faster by using the exponential map of spin
    # operators, expressed as a polynomial expansion,
    # http://www.emis.de/journals/SIGMA/2014/084/
    (evals, evecs) = eigen(n⋅S)
    return normalize(evecs[:, argmax(evals)])
end


struct System{N, N2}
    periodic::Bool # Periodic boundary conditions

    spin_rescaling::Float64 # Rescaling factor applied to all expectation values
    S::Vec3{SMatrix{N, N, ComplexF64, N2}} # Spin operators, [Sᵃ, Sᵇ] = i ϵᵃᵇᶜ Sᶜ
    
    L::Int # Chain length
    J::Float64 # Heisenberg coupling strength
    B::Vec3{Float64} # External field
    D::Float64 # Easy-axis strength, in the LL formalism
    Λ::Vector{SMatrix{N, N, ComplexF64, N2}} # Anisotropy matrices in GSD formalism

    rng::AbstractRNG
    ξ::Vector{SVector{N, ComplexF64}} # Random noise for Langevin integration 
    λ::Float64 # Damping coefficient

    Z::Vector{SVector{N, ComplexF64}} # Coherent state, taken to be normalized
    s::Vector{Vec3{Float64}} # Spin dipole expectation values

    ∇E::Vector{Vec3{Float64}} # Storage space for energy gradient, dE/d𝐬ᵢ

    # Temporary storage
    Z_temp1::Vector{SVector{N, ComplexF64}}
    Z_temp2::Vector{SVector{N, ComplexF64}}
    Z_temp3::Vector{SVector{N, ComplexF64}}
    Z_temp4::Vector{SVector{N, ComplexF64}}

    s_temp1::Vector{Vec3{Float64}}
    s_temp2::Vector{Vec3{Float64}}
    s_temp3::Vector{Vec3{Float64}}

    function System(; N, λ=0.1, periodic=false, spin_rescaling=1.0, L, J=0., B=[0., 0., 0.], D = 0.0, Λ=nothing, rng=nothing)

        S = spin_operators(N)

        if isnothing(rng)
            rng = MersenneTwister()
        end
        if isnothing(Λ)
            Λ = zeros(SMatrix{N, N, ComplexF64, N^2}, L)
        end
        if N > 0
            # Verify anisotropy matrices are Hermitian
            for i = 1:L
                @assert norm(Λ[i] - Λ'[i]) < 1e-12
            end
        end

        Z₀ = SVector{N,ComplexF64}((i==1 ? 1 : 0) for i=1:N)
        if N == 0
            s₀ = Vec3(0, 0, spin_rescaling)
        else
            s₀ = spin_rescaling * spin_bilinear(S, Z₀)
        end
        Z = fill(Z₀, L)
        ξ = fill(Z₀, L)
        s = fill(s₀, L)

        ∇E = zeros(Vec3{Float64}, L)

        Z_temp1 = zeros(SVector{N, ComplexF64}, L)
        Z_temp2 = zeros(SVector{N, ComplexF64}, L)
        Z_temp3 = zeros(SVector{N, ComplexF64}, L)
        Z_temp4 = zeros(SVector{N, ComplexF64}, L)

        s_temp1 = zeros(Vec3{Float64}, L)
        s_temp2 = zeros(Vec3{Float64}, L)
        s_temp3 = zeros(Vec3{Float64}, L)

        return new{N, N^2}(periodic, spin_rescaling, S, L, J, B, D, Λ, rng, ξ, λ, Z, s, ∇E, Z_temp1, Z_temp2, Z_temp3, Z_temp4, s_temp1, s_temp2, s_temp3)
    end
end

function normalize_spin(sys::System{N}, s) where N
    @assert N == 0
    return sys.spin_rescaling * normalize(s)
end

function set_state_from_dipole!(sys::System{N}, i, s::Vec3{T}) where {T, N}
    if N == 0
        sys.s[i] = normalize_spin(sys, s)
    else
        sys.Z[i] = infer_ket_from_dipole(sys.S, s)
        sys.s[i] = sys.spin_rescaling * spin_bilinear(sys.S, sys.Z[i])
    end
end

function set_expected_spins!(s, sys, Z)
    S = Ref(sys.S)
    @. s = sys.spin_rescaling * spin_bilinear(S, Z)
end

function local_energy_and_gradient(sys::System{N}, s, i) where {N}
    E = 0.
    ∇E = zero(Vec3{Float64})

    # Exchange interactions have the form J ∑_⟨i,j⟩ nᵢ nⱼ. Local energy
    # contribution at site i is then  (1/2) J nᵢ ∑ⱼ nⱼ, to handle double
    # counting of bonds.
    if i > 1 || sys.periodic
        i⁻ = mod1(i-1, sys.L)
        E += (1/2)*sys.J*(s[i⁻]⋅s[i])
        ∇E += sys.J*s[i⁻]
    end
    if i < sys.L || sys.periodic
        i⁺ = mod1(i+1, sys.L)
        E += (1/2)*sys.J*(s[i⁺]⋅s[i])
        ∇E += sys.J*s[i⁺]
    end

    # Zeeman coupling
    E += -sys.B ⋅ s[i]
    ∇E += -sys.B

    # LL easy axis, - D ∑ᵢ (nᵢᶻ)²
    E += sys.D*s[i][3]^2
    ∇E += Vec3(0, 0, 2sys.D*s[i][3])

    return (; E, ∇E)
end

function set_gradient!(∇E, sys, n)
    for i = 1:sys.L
        ∇E[i] = local_energy_and_gradient(sys, n, i).∇E
    end
end

function apply_local_hamiltonians!(out, sys, ∇E, Z)
    κ = sys.spin_rescaling
    for i = 1:length(Z)
        out[i] = κ * (∇E[i]⋅sys.S + sys.Λ[i]) * Z[i]
    end
    return nothing
end

function energy(sys::System{N}) where N
    E = 0.
    for i = 1:sys.L
        E += local_energy_and_gradient(sys, sys.s, i).E
        if N > 0
            E += sys.spin_rescaling * real(sys.Z[i]' * sys.Λ[i] * sys.Z[i])
        end
    end
    return E
end


function integrate_rhs_aux!(rhs, sys, Δt, kT, ∇E, Z)
    (; ξ, λ) = sys 
    apply_local_hamiltonians!(rhs, sys, ∇E, Z)
    for i = 1:sys.L
        rhs[i] = -im*√(2*Δt*kT*λ)*ξ[i] - Δt*(im+λ)*rhs[i]
    end
end

function integrate_rhs!(rhs, sys::System, Δt, kT, Z)
    set_expected_spins!(sys.s, sys, Z)
    set_gradient!(sys.∇E, sys, sys.s)
    integrate_rhs_aux!(rhs, sys, Δt, kT, sys.∇E, Z)
end

function heun_langevin_step!(sys::System, Δt, kT)
    rhs = sys.Z_temp1
    Z′ = sys.Z_temp2
    rhs′ = sys.Z_temp3

    randn!(sys.rng, sys.ξ)

    integrate_rhs!(rhs, sys, Δt, kT, sys.Z)
    @. Z′ = normalize(sys.Z + rhs)

    integrate_rhs!(rhs′, sys, Δt, kT, Z′)
    @. sys.Z = normalize(sys.Z + (rhs + rhs′)/2)

    set_expected_spins!(sys.s, sys, sys.Z)
end

function Random.rand!(sys::System)
    randn!(sys.rng, sys.Z)
    @. sys.Z /= norm(sys.Z)
    set_expected_spins!(sys.s, sys, sys.Z)
    nothing
end



## Code to generate figures and animation 

using Plots
using ColorSchemes
using LaTeXStrings
using JLD2
import Measures: mm
import Statistics: mean, std
import SpecialFunctions: erf

# pyplot() # Install PyPlot and uncomment for publication figures

include(joinpath("figures_1_and_2", "models_and_utils.jl"))
include(joinpath("figures_1_and_2", "fig1.jl"))
include(joinpath("figures_1_and_2", "fig2.jl"))
export fig1, fig2

include(joinpath("skyrmions", "LangevinSkyrmions.jl"))
fig3 = LangevinSkyrmions.fig3
quench_animation = LangevinSkyrmions.quench_animation
thermalize_animation = LangevinSkyrmions.thermalize_animation
export fig3, quench_animation, thermalize_animation


end # module


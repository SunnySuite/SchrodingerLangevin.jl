using GLMakie
using Sunny
import Statistics: mean
using LinearAlgebra
using JLD2
using Parameters
using Formatting
using DelimitedFiles

includet("model.jl")
includet("sim_utils.jl")
includet("plaquette_utils.jl")
includet("plot_utils.jl")



################################################################################
# Functions for processing Hao's data
################################################################################
raw_angles = readdlm("src/fig3/hao_angles.dat")

function angle_to_cs(angles)
    θ, ϕ, α1, α2 = angles
    return SA[
        exp(im*α1)*sin(θ)*cos(ϕ),
        cos(θ),
        exp(im*α2)*sin(θ)*sin(ϕ)
]
end

function reinterpret_hao(raw_angles)
    angles = reinterpret(SVector{4, Float64}, raw_angles)
    angles_new = zeros(SVector{4, Float64}, 10, 10, 1, 1)
    for (i, angle) ∈ enumerate(angles)
        r = fld1(i, 10)
        c = mod1(i, 10)
        angles_new[c, r, 1, 1] = angle # transposed 
    end
    map(angle_to_cs, angles_new)
end

function periodically_extend(gs, factor) 
    x, y = size(gs)[1:2]
    Zs = zeros(Sunny.CVec{3}, x*factor, y*factor, 1, 1)
    for i ∈ 1:factor, j ∈ 1:factor
        Zs[(i-1)*x+1:i*x, (j-1)*y+1:j*y, 1, 1] .= gs[:,:,1,1]
    end
    return Zs
end 

begin
    Zs = reinterpret_hao(raw_angles)
    Zs = periodically_extend(Zs, 10)
    dims_hao = (100, 100, 1)
    sys_hao = su3_skyrmion_model(dims_hao; D=19.0, h=15.349710144927537, rng)
    Sunny.init_from_coherents!(sys_hao, Zs)
    energy(sys_hao)/10000
end

# Generate and save trajectory
begin
    dims = (100,100,1)
    h = 15.35

    #= Set up system and initialize =#
    seed = 111
    rng = MersenneTwister(seed)
    sys = su3_skyrmion_model(dims; h, rng)

    raw_angles = readdlm("src/fig3/hao_angles.dat")
    Z₀ = reinterpret_hao(raw_angles)
    Z₀ = periodically_extend(Z₀, 10)
    println(size(Z₀))
    Sunny.init_from_coherents!(sys, Z₀)
end

# Run trial
begin
    kT = 0.11
    α = 0.1
    Δt = 0.01
    τs = [50.0, 100.0, 200.0, 300.0, 400.0]
    quench_dur = 100.0

    integrator = LangevinHeunP(sys, kT, α)
    steps_quench = round(Int, quench_dur/Δt)

    final_configs = []
    Es = []

    for τ ∈ τs
        println("Thermalization dur: τ")
        steps = round(Int, τ/Δt) 
        Sunny.init_from_coherents!(sys, Z₀)

        # Thermalize
        print("\tThermalizing...\t")
        integrator.kT = kT
        @time for _ ∈ 1:steps
            evolve!(integrator, Δt)
        end

        # Quench
        print("\tQuenching...   \t")
        integrator.kT = 0.0
        @time for _ ∈ 1:steps_quench
            evolve!(integrator, Δt)
        end

        # Save result and energy
        push!(final_configs, copy(sys._coherents))
        push!(Es, energy(sys))
    end

    Es ./= 10_000
end

# Plot results
begin
    vecs = sys.lattice.lat_vecs
    v₁, v₂ = vecs[:,1], vecs[:,2]
    χs = [plaquette_map(berry, Zs) for Zs ∈ final_configs]
    # χs = [round(Int, sum(abs.(χ))/2π) for χ ∈ χs]
    χs = [sum(abs.(χ))/2π for χ ∈ χs]
    texts = [format("τ = {:.1f}\nE = {:.4f} J₁\nχ = {:.1f}", τ, E, χ) for (τ, E, χ) in zip(τs, Es, χs)]
    fig, ax = plot_chirality_multi(final_configs, v₁, v₂;
        clims=(-1, 1), numcols = 1, texts
    )
    fig
end
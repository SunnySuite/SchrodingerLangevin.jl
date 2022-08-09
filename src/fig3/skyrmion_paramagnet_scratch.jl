using GLMakie
using Sunny
import Statistics: mean
using LinearAlgebra
using JLD2
using Parameters
using Formatting

includet("model.jl")
includet("sim_utils.jl")
includet("plaquette_utils.jl")
includet("plot_utils.jl")



function bin_idcs(Zs, idx; width = 10)
    num_steps = size(Zs)[end]
    i1 = max(1, idx - round(Int, width/2) + 1)
    i2 = min(num_steps, i1 + width)
    i1:i2
end

# Generate and save trajectory
begin
    step_dur = 10.0
    base = 1.15 
    kTs1 = [0.04*(base^k) for k ∈ 0:6]
    kTs2 = kTs1[end-1:-1:1] 
    kTs = vcat(kTs1, kTs2)
    push!(kTs, 0.0)
    dims = (100,100,1)
    h = 15.35
    Δt = 0.01
    α = 0.1

    #= Set up system and initialize =#
    seed = 111
    rng = MersenneTwister(seed)
    sys = su3_skyrmion_model(dims; h, rng)
    Z₀ = fill(SVector{3, ComplexF64}([0, 1, 0]), size(sys._coherents))
    Sunny.init_from_coherents!(sys, Z₀)



    num_steps = round(Int, step_dur/Δt)
    integrator = LangevinHeunP(sys, kTs[1], α)

    trajectory = zeros(SVector{3, ComplexF64}, size(sys._coherents)..., num_steps*length(kTs))
    kTs_all = zeros(num_steps*length(kTs))

    c = 1
    for (i, kT) ∈ enumerate(kTs)
        println(kT)
        integrator.kT = kT
        @time for _ ∈ 1:num_steps
            evolve!(integrator, Δt)
            trajectory[:,:,:,:,c] .= sys._coherents
            kTs_all[c] = kT
            c += 1
        end
    end
end

begin
    χs = zeros(size(trajectory)[end])
    for i ∈ 1:size(trajectory)[end]
        χ = plaquette_map(berry, trajectory[:,:,:,:,i])
        χs[i] = round(Int, sum(χ)/2π)
    end
end


begin
    text = ["Total charge: $(χ)\nT=$T" for (χ, T) ∈ zip(χs, kTs_all)]
    lat_vecs = sys.lattice.lat_vecs
    v1, v2 = lat_vecs[:,1], lat_vecs[:,2]
    animate_chirality(trajectory, v1, v2;
        filename="nucleation.mp4", clims=(-1, 1), skip_interval=20,
        text
    )
end

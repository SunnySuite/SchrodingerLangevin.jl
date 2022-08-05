using GLMakie
using Sunny
import Statistics: mean
using LinearAlgebra
using JLD2

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

dur = 200.0
dims = (100,100,1)
h = 15.35
Δt = 0.01
α = 0.1

#= Set up system and initialize =#
seed = 111
rng = MersenneTwister(seed)
sys = su3_skyrmion_model(dims; h, rng)
rand!(sys)

#= Run and save trajectory =#
integrator = LangevinHeunP(sys, kT, α)
Zs = ket_trajectory!(integrator, Δt, dur)

#= Save results =#
save("trajectory.jld2", "Zs", Zs, "Δt", Δt, "h", h)


τ = 0.1
base = 10.0
τs = [τ*(base^k) for k \include     
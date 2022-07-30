include("SchrodingerLangevin.jl")
using .SchrodingerLangevin

using Random
using Plots

rng = MersenneTwister(1)
N = 2
L = 2
kT = 0.0
sys = SchrodingerLangevin.System(; N, spin_rescaling=1.0, J=1.0, L, rng)

Δt = 0.01
randn!(sys.Z)

dur = 100.0
n_steps = round(Int, dur/Δt)
Es = zeros(n_steps)

for i ∈ 1:n_steps
    SchrodingerLangevin.heun_langevin_step!(sys, Δt, kT)
    Es[i] = SchrodingerLangevin.energy(sys)
end

plot(Es)
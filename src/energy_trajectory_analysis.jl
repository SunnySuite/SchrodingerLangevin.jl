includet("SchrodingerLangevin.jl")
using .SchrodingerLangevin

using Random
# using Plots
using LinearAlgebra
using FFTW
using GLMakie 
using Interpolations
using Statistics

rng = MersenneTwister(111)

function dimer_exchange(J)
    v = [0.0 1.0 -1.0 0.0] / √2
    return (J/4)*I(4) - J*v'*v
end

function energy_trajectory!(sys, dur, Δt, kT)
    n = round(Int, dur/Δt)
    Es = zeros(n+1)
    Es[1] = energy(sys)
    for i ∈ 1:n
        heun_langevin_step!(sys, Δt, kT)
        Es[i+1] = energy(sys)
    end
    return Es
end

function classical_system(; J=1.0, rng=nothing)
    System(; N=2, L=2, J, rng)
end

function entangled_pair(; J=1.0, rng=nothing)
    Λ = dimer_exchange(J)
    System(; N=4, L=1, J=0.0, Λ=[Λ], rng)
end

function generate_statistics(Δt, num_samples, kTs;
    sys_func, 
    bin_func,
    dur_burnin=100.0,
)
    μs = zero(kTs)
    σs = zero(kTs)
    for (i, kT) in enumerate(kTs)
        println("Collecting statistics for kT=$kT...")

        sys = sys_func()
        rand!(sys)
        energy_trajectory!(sys, dur_burnin, Δt, kT)
        dur = bin_func(kT)

        Es = zeros(num_samples)
        for j ∈ 1:num_samples
            E = energy_trajectory!(sys, dur, Δt, kT)
            Es[j] = mean(E) 
        end
        μs[i] = mean(Es)
        σs[i] = std(Es)
    end
    return (; μs, σs)
end

function energy_classical(β, J=1.0)
    (1/β) - (J/4)*coth(β*J/4)
end

function energy_entangled(β, J=1.0)
    if β >= 1000 || isinf(β)
        return J > 0 ? -0.75 : -0.25
    end
    num = 24 + 18J*β + 6(J*β)^2 + (J*β)^3 + 6*exp(J*β)*(J*β - 4)
    denom = 4*β*(2 - 2*exp(J*β) + 2*J*β + (J*β)^2) 
    return num/denom
end

function energy_quantum(β, J=1.0)
    if β >= 1000 || isinf(β)
        return J > 0 ? -0.75 : -0.25
    end
    -3J*(exp(β*J) - 1) / (4*(exp(β*J) + 3))
end

# Parameters
begin
    J = -1.0
    Δt = 0.01
    dur = 500.0
    dur_bi = 100.0
    kTs = [0.01, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 2.5, 2.75, 3.0]
    num_samples = 250
end

begin
    Es_all = []
    for kT ∈ kTs
        sys_e = entangled_pair(; J)
        rand!(sys_e)
        Es_local = []
        @time for i ∈ 1:num_samples
            rand!(sys_e)
            energy_trajectory!(sys_e, dur_bi, Δt, kT) # Burnin
            Es = energy_trajectory!(sys_e, dur, Δt, kT) # Collect data
            push!(Es_local, Es)
        end
        push!(Es_all, Es_local)
    end
end

begin
    ts = collect(0:length(Es_all[1][1])-1) .* Δt
    i_max = findfirst(x -> x > 50, ts)
    fig = GLMakie.Figure(resolution=(800,600))
    num_cols = 3
    for (i, Es) ∈ enumerate(Es_all)
        ax = GLMakie.Axis(fig[fldmod1(i, num_cols)...];
            title = "kT=$(kTs[i])",
        )
        es = fft.(Es)
        acs = [real(ifft(e .* conj.(e))) for e ∈ es]
        ac = sum(acs) ./ length(acs)
        lines!(ax, ts[1:i_max], ac[1:i_max])
    end
    fig
end

# begin
#     J = 1.0
#     kTs = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
#     # kTs = 0.1:0.1:1.0 
#     bin_dur = [40.0, 25.0, 10.0, 6.0, 5.0, 5.0, 5.0]
#     itp = LinearInterpolation(kTs, bin_dur)
#     bin_func = kT -> kT < 0.1 ? 1.0 : itp(kT)
# end


# Entangled AFM (J=1.0)
kTs =      [0.01, 0.25, 0.5,  0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 2.5, 2.75, 3.0]
bin_durs = [40.0, 35.0, 20.0, 10.0, 9.0, 7.5,  5.0, 5.0,  5.0, 4.0, 3.0, 3.0,  3.0]
itp = LinearInterpolation(kTs, bin_durs)
bin_en_afm = kT -> kT < 0.01 ? 1.0 : itp(kT)

# Entangled FM (J=-1.0)
kTs =      [0.01, 0.25, 0.5,  0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 2.5, 2.75, 3.0]
bin_durs = [40.0, 25.0, 15.0, 10.0, 7.5, 5.0,  5.0, 5.0,  4.0, 3.0, 3.0, 3.0,  3.0]
itp = LinearInterpolation(kTs, bin_durs)
bin_en_fm = kT -> kT < 0.01 ? 1.0 : itp(kT)

# Classical AFM (J=1.0)
kTs =      [0.01, 0.25, 0.5,  0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 2.5, 2.75, 3.0]
bin_durs = []
itp = LinearInterpolation(kTs, bin_durs)
bin_en_afm = kT -> kT < 0.01 ? 1.0 : itp(kT)

# Classical FM (J=-1.0)
kTs =      [0.01, 0.25, 0.5,  0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 2.5, 2.75, 3.0]
bin_durs = []
itp = LinearInterpolation(kTs, bin_durs)
bin_en_afm = kT -> kT < 0.01 ? 1.0 : itp(kT)
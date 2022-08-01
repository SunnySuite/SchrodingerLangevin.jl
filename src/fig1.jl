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


function fig1()
    #= Set sampling bin sizes =#
    kTs = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    bin_dur = [40.0, 25.0, 10.0, 6.0, 5.0, 5.0, 5.0]
    itp = LinearInterpolation(kTs, bin_dur)
    bin_func = kT -> kT < 0.1 ? 1.0 : itp(kT)

    #= Trial parameters =#
    Δt = 0.1
    num_samples = 1000
    kTs = 0.0:0.1:2.0

    #= Generate data for classical and entangled unit models =#
    rng = MersenneTwister(111)

    sys_func = () -> entangled_pair(; J=1.0, rng) 
    (; μs, σs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func)
    μs_afm_en, σs_afm_en = μs, σs

    sys_func = () -> entangled_pair(; J=-1.0, rng) 
    (; μs, σs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func)
    μs_fm_en, σs_fm_en = μs, σs

    sys_func = () -> classical_system(; J=1.0, rng) 
    (; μs, σs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func)
    μs_afm_cl, σs_afm_cl = μs, σs

    sys_func = () -> classical_system(; J=-1.0, rng) 
    (; μs, σs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func)
    μs_fm_cl, σs_fm_cl = μs, σs


    #= Plot results =#
    kTs_ref = 0.000:0.001:2.0

    # FM
    J = -1.0
    E_fm_cl = map(kT -> energy_classical(1/kT, J), kTs_ref)
    E_fm_qu = map(kT -> energy_quantum(1/kT, J), kTs_ref)
    E_fm_en = map(kT -> energy_entangled(1/kT, J), kTs_ref) 
    p1 = plot(kTs_ref, E_fm_cl;
        label="LL",
        linestyle=:dash,
        color = :black,
        ylabel = "E",
    )
    plot!(kTs_ref, E_fm_en;
        label="Entangled Unit",
        linestyle=:dashdot,
        color = :black,
    )
    plot!(kTs_ref, E_fm_qu;
        label="Quantum",
        linestyle=:dashdotdot,
        color = :black,
    )
    scatter!(kTs, μs_fm_cl;
        label = "LL",
        markershape = :circle,
        color = 1,
    )
    scatter!(kTs, μs_fm_en;
        markershape = :diamond,
        label = "Entangled Unit",
        color = 2,
    )

    ## AFM
    J = 1.0
    E_afm_cl = map(kT -> energy_classical(1/kT, J), kTs_ref)
    E_afm_qu = map(kT -> energy_quantum(1/kT, J), kTs_ref)
    E_afm_en = map(kT -> energy_entangled(1/kT, J), kTs_ref) 
    p2 = plot(kTs_ref, E_afm_cl;
        label="LL",
        linestyle = :dash,
        color = :black,
        legend=false,
        xlabel = "kT",
        ylabel = "E",
    )
    plot!(kTs_ref, E_afm_en; 
        label="Entangled Unit",
        linestyle=:dashdot,
        color = :black,
    )
    plot!(kTs_ref, E_afm_qu;
        label="Quantum",
        linestyle=:dashdotdot,
        color = :black,
    )
    scatter!(kTs, μs_afm_cl;
        color = 1,
        markershape = :circle,
    )
    scatter!(kTs, μs_afm_en;
        markershape = :diamond,
        color = 2,
    )

    layout = @layout [a; b]
    return plot(p1, p2;
        layout,
        title = ["FM" "AFM"],
    )
end
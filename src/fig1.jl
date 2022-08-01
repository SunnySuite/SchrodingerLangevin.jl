# using Plots
# using ColorSchemes
# using LaTeXStrings
# using Interpolations
# import Measures: mm
# import Statistics: mean, std
# using JLD2


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

function classical_pair(; J=1.0, rng=nothing)
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
    #= Set sampling bin widths (decorrelation times) for different models and temperatures =#
    # Entangled AFM (J=1.0)
    kTs =      [0.01, 0.25, 0.5,  0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]

    bin_durs = [40.0, 35.0, 20.0, 10.0, 9.0, 7.5,  5.0, 5.0,  5.0, 4.0, 3.0, 3.0,  3.0]
    itp = LinearInterpolation(kTs, bin_durs)
    bin_en_afm = kT -> kT < 0.01 ? 1.0 : itp(kT)

    # Entangled FM (J=-1.0)
    bin_durs = [40.0, 25.0, 15.0, 10.0, 7.5, 5.0,  5.0, 5.0,  4.0, 3.0, 3.0, 3.0,  3.0]
    itp = LinearInterpolation(kTs, bin_durs)
    bin_en_fm = kT -> kT < 0.01 ? 1.0 : itp(kT)

    # Classical AFM (J=1.0)
    bin_durs = [40.0, 25.0, 12.0, 8.5,  8.0, 7.5,  5.0, 4.0,  4.0, 3.0, 3.0, 3.0,  2.5]
    itp = LinearInterpolation(kTs, bin_durs)
    bin_cl_afm = kT -> kT < 0.01 ? 1.0 : itp(kT)

    # Classical FM (J=-1.0)
    bin_durs = [40.0, 30.0, 15.0, 12.0, 10.0, 7.5,  5.0, 4.9,  4.0, 3.5, 2.75, 2.5,  2.45]
    itp = LinearInterpolation(kTs, bin_durs)
    bin_cl_fm = kT -> kT < 0.01 ? 1.0 : itp(kT)


    #= Trial parameters =#
    kT_max = 1.5
    Δt = 0.1
    num_samples = 10000
    kTs = 0.0:0.1:kT_max


    #= Generate data for classical and entangled unit models =#
    rng = MersenneTwister(111)

    println("Collecting statistics to pair of SU(2) coherent spins (J=1)")
    sys_func = () -> entangled_pair(; J=1.0, rng) 
    (; μs, σs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func=bin_en_afm)
    μs_afm_en, σs_afm_en = μs, σs

    println("Collecting statistics for single SU(4) engtangled pair (J=1)") 
    sys_func = () -> entangled_pair(; J=-1.0, rng) 
    (; μs, σs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func=bin_en_fm)
    μs_fm_en, σs_fm_en = μs, σs

    println("Collecting statistics to pair of SU(2) coherent spins (J=-1)")
    sys_func = () -> classical_pair(; J=1.0, rng) 
    (; μs, σs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func=bin_cl_afm)
    μs_afm_cl, σs_afm_cl = μs, σs

    println("Collecting statistics for single SU(4) engtangled pair (J=-1)") 
    sys_func = () -> classical_pair(; J=-1.0, rng) 
    (; μs, σs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func=bin_cl_fm)
    μs_fm_cl, σs_fm_cl = μs, σs


    #= Calculate analytical energies =#
    kTs_ref = 0.000:0.001:kT_max

    J = -1.0
    E_fm_cl = map(kT -> energy_classical(1/kT, J), kTs_ref)
    E_fm_qu = map(kT -> energy_quantum(1/kT, J), kTs_ref)
    E_fm_en = map(kT -> energy_entangled(1/kT, J), kTs_ref) 

    J = 1.0
    E_afm_cl = map(kT -> energy_classical(1/kT, J), kTs_ref)
    E_afm_qu = map(kT -> energy_quantum(1/kT, J), kTs_ref)
    E_afm_en = map(kT -> energy_entangled(1/kT, J), kTs_ref) 


    #= Plot results =#
    params = (;
        xlabelfontsize = 16,
        ylabelfontsize = 14,
        legendfontsize = 14,
        xtickfontsize = 12,
        ytickfontsize = 12,
        palette=:seaborn_colorblind,
    )
    marker1 = (;
        color = 1,
        markershape = :cross,
        markersize = 5.0,
        markerstrokewidth = 2.0,
    )
    marker2 = (;
        color = 2,
        markershape = :xcross,
        markersize = 4.0,
        markerstrokewidth = 2.0,
    )
    line_cl = (;
        linewidth=1.5,
        linestyle = :dash,
        color = :black,
    )
    line_en = (;
        linestyle=:dot,
        linewidth=2.0,
        color = :black,
    )
    line_qu = (;
        linewidth=1.5,
        color = :black,
    )
    yticks_fm  = 0.0:-0.05:-0.25
    yticks_afm = 0.0:-0.125:-0.75
    xticks = 0.0:0.5:kTs[end]

    p1 = plot(; 
        yticks = (yticks_fm, [L"%$y" for y ∈ yticks_fm]),
        ylims = (-0.25*1.05, 0.0),
        xticks = (xticks, [L"%$x" for x ∈ xticks]),
        params...
    )
    plot!(kTs_ref, E_fm_cl;
        label="SU(2)",
        ylabel = L"E",
        legend=false,
        line_cl...
    )
    plot!(kTs_ref, E_fm_en;
        label="SU(4)",
        color = :black,
        line_en...
    )
    plot!(kTs_ref, E_fm_qu;
        label="Quantum",
        color = :black,
        line_qu...
    )
    scatter!(kTs, μs_fm_cl;
        label = "LL",
        marker1...
    )
    scatter!(kTs, μs_fm_en;
        label = "Entangled Unit",
        marker2...
    )

    ## AFM
    p2 = plot(;
        yticks = (yticks_afm, [L"%$y" for y ∈ yticks_afm]),
        ylims = (-0.75*1.05, 0.0),
        xticks = (xticks, [L"%$x" for x ∈ xticks]),
        legend=:bottomright,
        params...
    )
    plot!(kTs_ref, E_afm_cl;
        label="SU(2)",
        color = :black,
        xlabel = L"k_bT",
        ylabel = L"E",
        line_cl...
    )
    plot!(kTs_ref, E_afm_en; 
        label="SU(4)",
        line_en...
    )
    plot!(kTs_ref, E_afm_qu;
        label="Quantum",
        linewidth=1.0,
        color = :black,
        line_qu...
    )
    scatter!(kTs, μs_afm_cl;
        label="SU(2) numerical",
        marker1...
    )
    scatter!(kTs, μs_afm_en;
        label="SU(4) numerical",
        marker2...
    )

    layout = @layout [a; b]
    p = plot(p1, p2;
        layout,
        title = [L"J=-1.0" L"J=1.0"],
    )
    plot!(;size=(600,600))
    return p
end
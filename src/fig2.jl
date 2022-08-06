## Just for development. Delete these comments before publishing.
include("SchrodingerLangevin.jl")
include("models_and_utils.jl")
using .SchrodingerLangevin
using Plots
using ColorSchemes
using LaTeXStrings
using Interpolations
using LinearAlgebra
using Random
import Measures: mm
import Statistics: mean, std


begin 
function fig2()
    #= Generate data for classical and entangled unit models =#

    # Trial parameters
    kT_max = 1.5 # Maximum temperature
    kTs = 0.0:0.1:kT_max
    Δt = 0.01  # Integration time step
    α = 1.0  # Empirical damping parameter
    num_samples = 10  # Number of estimates of the energy (number of trajectories)
    num_samples_fm = 100  # Number of estimates of the energy for J=-1.0
    dur_burnin = 10.0  # Burnin-duration for each trajectory
    dur_trajectory = 20.0  # Duration of each trajectory
    rng = MersenneTwister(17)

    # Run simulations
    println("\nCollecting statistics to pair of SU(2) coherent spins (J=1)")
    sys_func = () -> entangled_pair(; J=1.0, rng, α) 
    (; μs, sems) = generate_statistics(Δt, num_samples, kTs; sys_func, dur_burnin, dur_trajectory)
    μs_afm_en, sems_afm_en = μs, sems

    println("\nCollecting statistics for single SU(4) engtangled pair (J=1)") 
    sys_func = () -> entangled_pair(; J=-1.0, rng, α) 
    (; μs, sems) = generate_statistics(Δt, num_samples_fm, kTs; sys_func, dur_burnin, dur_trajectory)
    μs_fm_en, sems_fm_en = μs, sems

    println("\nCollecting statistics to pair of SU(2) coherent spins (J=-1)")
    sys_func = () -> classical_pair(; J=1.0, rng, α) 
    (; μs, sems) = generate_statistics(Δt, num_samples, kTs; sys_func, dur_burnin, dur_trajectory)
    μs_afm_cl, sems_afm_cl = μs, sems

    println("\nCollecting statistics for single SU(4) engtangled pair (J=-1)") 
    sys_func = () -> classical_pair(; J=-1.0, rng) 
    (; μs, sems) = generate_statistics(Δt, num_samples_fm, kTs; sys_func, dur_burnin, dur_trajectory)
    μs_fm_cl, sems_fm_cl = μs, sems


    #= Calculate analytical energies =#
    kTs_ref = 0.000:0.001:kT_max
    J = -1.0
    E_fm_cl = map(kT -> energy_pair_dipole(1/kT, J), kTs_ref)
    E_fm_en = map(kT -> energy_pair_entangled(1/kT, J), kTs_ref) 
    E_fm_qu = map(kT -> energy_pair_quantum(1/kT, J), kTs_ref)
    J = 1.0
    E_afm_cl = map(kT -> energy_pair_dipole(1/kT, J), kTs_ref)
    E_afm_en = map(kT -> energy_pair_entangled(1/kT, J), kTs_ref) 
    E_afm_qu = map(kT -> energy_pair_quantum(1/kT, J), kTs_ref)


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
        label = "Dipole numerical",
        color = 1,
        alpha = 0.65,
        markersize = 7.0,
        markerstrokewidth = 0.5,
    )
    marker2 = ( ;
        label = "SU(4) numerical",
        color = 2,
        alpha = 0.65,
        markersize = 7.0,
        markerstrokewidth = 0.5,
    )
    error_bars = (;
        markersize=8.0,
        markerstrokewidth=1.0,
        markerstrokecolor=:black,
        markeralpha=1.00,
        linewidth = 0.0,
        label=false
    )
    line_cl = (;
        label = "Dipole",
        linewidth = 1.0,
        linestyle = :dash,
        linealpha=0.6,
        color = :black,
    )
    line_en = (;
        label = "SU(4)",
        linestyle = :dot,
        linewidth = 1.5,
        linealpha=0.6,
        color = :black,
    )
    line_qu = (;
        label = "Quantum",
        linewidth=1.0,
        linealpha=0.6,
        color = :black,
    )
    yticks_fm  = 0.0:-0.05:-0.25
    yticks_afm = 0.0:-0.125:-0.75
    xticks = 0.0:0.5:kTs[end]

    # FM (J < 0)
    p1 = plot(; 
        yticks = (yticks_fm, [L"%$y" for y ∈ yticks_fm]),
        ylims = (-0.25*1.05, 0.05/3),
        xticks = (xticks, [L"%$x" for x ∈ xticks]),
        legend=false,
        params...
    )
    plot!(kTs_ref, E_fm_cl; ylabel = L"$\left\langle E \right\rangle$", line_cl...)
    plot!(kTs_ref, E_fm_en; line_en...)
    plot!(kTs_ref, E_fm_qu; line_qu...)
    plot!(kTs, μs_fm_cl; yerr=sems_fm_cl, error_bars...)
    scatter!(kTs, μs_fm_cl; marker1...)
    plot!(kTs, μs_fm_en; yerr=sems_fm_en, error_bars...)
    scatter!(kTs, μs_fm_en; marker2...)

    # AFM (J > 0)
    p2 = plot(;
        yticks = (yticks_afm, [L"%$y" for y ∈ yticks_afm]),
        ylims = (-0.75*1.05, 0.05),
        xticks = (xticks, [L"%$x" for x ∈ xticks]),
        legend=:bottomright,
        params...
    )
    plot!(kTs_ref, E_afm_cl; xlabel = L"k_bT", ylabel = L"$\left\langle E \right\rangle$", line_cl...)
    plot!(kTs_ref, E_afm_en; line_en...)
    plot!(kTs_ref, E_afm_qu; line_qu...)
    plot!(kTs, μs_afm_cl; yerr=sems_afm_cl, error_bars...)
    scatter!(kTs, μs_afm_cl; marker1...)
    plot!(kTs, μs_afm_en; yerr=sems_afm_en, error_bars...)
    scatter!(kTs, μs_afm_en; marker2...)

    layout = @layout [a; b]
    p = plot(p1, p2;
        layout,
        size=(600,700),
        title = [L"(\mathrm{a})\quad J=-1.0" L"(\mathrm{b})\quad J=1.0"],
        left_margin=3mm,
    )

    savefig("fig2.pdf")

    return p
end
fig2()
end
includet("SchrodingerLangevin.jl")
includet("models_and_utils.jl")
using .SchrodingerLangevin
using Plots
pyplot()
using ColorSchemes
using LaTeXStrings
using Interpolations
using Random
import Measures: mm
import Statistics: mean, std

function fig2()
    #= Generate data for classical and entangled unit models =#
    # Trial parameters
    kT_max = 2.0
    kTs = 0.0:0.15:kT_max
    Δt = 0.1
    num_samples = 10000
    bin_func = kT -> 5.0 # Sampling bin width. Uniform value sufficient for estimating mean.
    rng = MersenneTwister(111)

    D = -1.0
    println("\nCollecting statistics for SU(2) spin with traditional anisotropy.")
    sys_func = () -> dipole_anisotropy(; D, rng) 
    (; μs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func)
    μs_cl = μs

    println("\nCollecting statistics for SU(3) spin with SU(3) anisotropy.")
    sys_func = () -> su3_anisotropy(; D, rng) 
    (; μs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func)
    μs_su3 = μs



    #= Calculate analytical energies =#
    kTs_ref = 0.000:0.001:kT_max
    E_cl = map(kT -> energy_aniso_dipole(1/kT, D), kTs_ref)
    E_su3 = map(kT -> energy_aniso_su3(1/kT, D), kTs_ref) 
    E_qu = map(kT -> energy_aniso_quantum(1/kT, D), kTs_ref)


    #= Plot results =#
    params = (;
        xlabelfontsize = 16,
        ylabelfontsize = 14,
        legendfontsize = 14,
        xtickfontsize = 12,
        ytickfontsize = 12,
        palette=:seaborn_colorblind,
    )
    marker1 = ( ;
        label = "Dipole aniso numerical",
        color = 1,
        markershape = :cross,
        markersize = 7.0,
        markerstrokewidth = 4.0,
    )
    marker2 = (;
        label = "SU(3) aniso numerical",
        color = 2,
        markershape = :xcross,
        markersize = 5.0,
        markerstrokewidth = 3.7,
    )
    line_cl = (;
        label = "Dipole aniso",
        linewidth = 1.5,
        linestyle = :dash,
        color = :black,
    )
    line_su3 = (;
        label = "SU(3) aniso",
        linestyle = :dot,
        linewidth = 2.0,
        color = :black,
    )
    line_qu = (;
        label = "Quantum",
        linewidth=1.5,
        color = :black,
    )
    yticks = -0.4:-0.1:-1.20
    xticks = 0.0:0.5:kTs[end]

    # FM (J < 0)
    p = plot(; 
        yticks = (yticks, [L"%$y" for y ∈ yticks]),
        xticks = (xticks, [L"%$x" for x ∈ xticks]),
        ylims = (-1.30, -0.35),
        xlims = (-0.05, 2.05),
        xlabel = L"E",
        ylabel = L"k_bT",
        legend=:bottomright,
        size=(600,350),
        params...
    )
    plot!(kTs_ref, E_cl; line_cl...)
    plot!(kTs_ref, E_su3; line_su3...)
    plot!(kTs_ref, E_qu; line_qu...)
    scatter!(kTs, μs_cl; marker1...)
    scatter!(kTs, μs_su3; marker2...)

    return p
end
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
    kT_max = 1.5
    kTs = 0.0:0.1:kT_max
    Δt = 0.1
    num_samples = 20_000
    bin_func = kt -> 10.0 # Sampling bin width. Uniform value sufficient for estimating mean.
    rng = MersenneTwister(111)

    # Run simulations
    println("\nCollecting statistics to pair of SU(2) coherent spins (J=1)")
    sys_func = () -> entangled_pair(; J=1.0, rng) 
    (; μs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func)
    μs_afm_en = μs

    println("\nCollecting statistics for single SU(4) engtangled pair (J=1)") 
    sys_func = () -> entangled_pair(; J=-1.0, rng) 
    (; μs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func)
    μs_fm_en = μs

    println("\nCollecting statistics to pair of SU(2) coherent spins (J=-1)")
    sys_func = () -> classical_pair(; J=1.0, rng) 
    (; μs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func)
    μs_afm_cl = μs

    println("\nCollecting statistics for single SU(4) engtangled pair (J=-1)") 
    sys_func = () -> classical_pair(; J=-1.0, rng) 
    (; μs) = generate_statistics(Δt, num_samples, kTs; sys_func, bin_func)
    μs_fm_cl = μs


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
    marker2 = ( ;
        label = "SU(4) numerical",
        color = 2,
        alpha = 0.65,
        # markershape = :cross,
        markersize = 7.0,
        markerstrokewidth = 1.0,
    )
    marker1 = (;
        label = "Dipole numerical",
        color = 1,
        alpha = 0.65,
        # markershape = :xcross,
        markersize = 7.0,
        markerstrokewidth = 1.0,
    )
    line_cl = (;
        label = "Dipole",
        linewidth = 2.0,
        linestyle = :dash,
        color = :black,
    )
    line_en = (;
        label = "SU(4)",
        linestyle = :dot,
        linewidth = 2.5,
        color = :black,
    )
    line_qu = (;
        label = "Quantum",
        linewidth=2.0,
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
    scatter!(kTs, μs_fm_cl; marker1...)
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
    scatter!(kTs, μs_afm_cl; marker1...)
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
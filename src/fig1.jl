## Just for development. Delete these comments before publishing.
# using Revise
include("SchrodingerLangevin.jl")
include("models_and_utils.jl")
using .SchrodingerLangevin
using Plots
using ColorSchemes
using LaTeXStrings
using Interpolations
using Random
import Measures: mm
import Statistics: mean
using Measurements

pyplot()

begin 
function fig1()
    #= Generate data for classical and entangled unit models =#
    # Trial parameters
    kT_max = 2.0  # Maximum temperature
    kTs = range(0.0, kT_max, length=16)
    Δt = 0.01  # Integration time step
    α = 1.0  # Empirical damping parameter
    num_samples = 100  # Number of estimates of the energy (number of trajectories)
    dur_burnin = 10.0  # Burnin-duration for each trajectory
    dur_trajectory = 10.0  # Duration of each trajectory
    D = -1.0  # Onsite anisotropy
    rng = MersenneTwister(12)

    println("\nCollecting statistics for SU(2) spin with traditional anisotropy.")
    sys_func = () -> dipole_anisotropy(; D, rng, α) 
    (; μs, sems) = generate_statistics(Δt, num_samples, kTs; sys_func, dur_trajectory, dur_burnin)
    μs_cl, sems_cl = μs, sems

    println("\nCollecting statistics for SU(3) spin with SU(3) anisotropy.")
    sys_func = () -> su3_anisotropy(; D, rng, α) 
    (; μs, sems) = generate_statistics(Δt, num_samples, kTs; sys_func, dur_trajectory, dur_burnin)
    μs_su3, sems_su3 = μs, sems



    #= Calculate analytical energies =#
    kTs_ref = 0.000:0.001:kT_max
    E_cl = map(kT -> energy_aniso_dipole(1/kT, D), kTs_ref)
    E_su3 = map(kT -> energy_aniso_su3(1/kT, D), kTs_ref) 
    E_qu = map(kT -> energy_aniso_quantum(1/kT, D), kTs_ref)


    #= Plot results =#
    yticks = -0.4:-0.1:-1.00
    xticks = 0.0:0.5:kT_max
    params = (;
        xlabelfontsize = 16,
        ylabelfontsize = 14,
        legendfontsize = 14,
        xtickfontsize = 12,
        ytickfontsize = 12,
        palette=:RdBu,
        yticks = (yticks, [L"%$y" for y ∈ yticks]),
        xticks = (xticks, [L"%$x" for x ∈ xticks]),
        ylims = (-1.05, -0.35),
        xlims = (-0.05, 2.05),
        ylabel = L"$\left\langle E \right\rangle$",
        xlabel = L"$k_b T$",
        legend=:bottomright,
        size=(600,425),
        left_margin=3mm,
        right_margin=3mm,
    )
    error_bars = (;
        markersize=10.0,
        markerstrokewidth=1.75,
        markeralpha=0.9,
        linewidth = 0.0,
        label=false
    )
    line_cl = (;
        label = "Classical (dipole only)",
        linewidth = 1.25,
        linestyle = :dash,
        linealpha=0.6,
        color = :black,
    )
    line_su3 = (;
        label = "Classical (SU(3) spin)",
        linestyle = :dot,
        linewidth = 1.25,
        linealpha=0.6,
        color = :black,
    )
    line_qu = (;
        label = "Quantum",
        linewidth=1.25,
        linealpha=0.6,
        color = :black,
    )

    p = plot(; params...)

    plot!(kTs_ref, E_cl; line_cl...)
    plot!(kTs_ref, E_su3; line_su3...)
    plot!(kTs_ref, E_qu; line_qu...)

    plot!(kTs, μs_su3;
        yerr=sems_su3,
        markerstrokecolor=10,
        error_bars...
    )
    plot!(kTs, μs_cl;
       yerr=sems_cl,
       markerstrokecolor=3,
       error_bars...
    )

    savefig("fig1.pdf")

    return p
end

fig1()
end
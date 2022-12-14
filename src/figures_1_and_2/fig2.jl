function fig2()
    #= Generate data for classical and entangled unit models =#

    # Trial parameters
    kT_max = 1.5 # Maximum temperature
    kTs = 0.0:0.1:kT_max
    Δt = 0.01  # Integration time step
    λ = 1.0  # Empirical damping parameter
    num_samples = 10  # Number of estimates of the energy (number of trajectories)
    dur_burnin = 10.0  # Burnin-duration for each trajectory
    dur_trajectory = 100.0  # Duration of each trajectory
    rng = MersenneTwister(12)

    # Run simulations
    println("\nCollecting statistics to pair of SU(2) coherent spins (J=-1)")
    sys_func = () -> classical_pair(; J=-1.0, rng, λ) 
    (; μs, sems) = generate_statistics(Δt, num_samples, kTs; sys_func, dur_burnin, dur_trajectory)
    μs_fm_cl, sems_fm_cl = μs, sems

    println("\nCollecting statistics to pair of SU(4) entangled pair (J=-1)")
    sys_func = () -> entangled_pair(; J=-1.0, rng, λ) 
    (; μs, sems) = generate_statistics(Δt, num_samples, kTs; sys_func, dur_burnin, dur_trajectory)
    μs_fm_en, sems_fm_en = μs, sems

    println("\nCollecting statistics for single SU(2) coherent spins (J=1)") 
    sys_func = () -> classical_pair(; J=1.0, rng, λ) 
    (; μs, sems) = generate_statistics(Δt, num_samples, kTs; sys_func, dur_burnin, dur_trajectory)
    μs_afm_cl, sems_afm_cl = μs, sems

    println("\nCollecting statistics for single SU(4) engtangled pair (J=1)") 
    sys_func = () -> entangled_pair(; J=1.0, rng, λ) 
    (; μs, sems) = generate_statistics(Δt, num_samples, kTs; sys_func, dur_burnin, dur_trajectory)
    μs_afm_en, sems_afm_en = μs, sems


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
    dipole_color = :blue
    su4_color = :black
    params = (;
        xlabelfontsize = 16,
        ylabelfontsize = 14,
        legendfontsize = 14,
        xtickfontsize = 12,
        ytickfontsize = 12,
        palette=:RdBu,
    )
    error_bars = (;
        markersize=10.0,
        markerstrokewidth=1.0,
        markeralpha=0.6,
        linewidth = 0.0,
        label=false,
    )
    line_cl = (;
        label = false, 
        linewidth = 1.5,
        linestyle = :dot,
        linealpha=0.9,
        color = dipole_color,
    )
    line_cl_label = (;
        label = "Classical (dipoles only)",
        linewidth = 1.5,
        linestyle = :dot,
        linealpha=1.0,
        color = dipole_color,
    )
    line_en = (;
        label = "Classical (SU(4) spin)",
        linestyle = :dash,
        linewidth = 1.0,
        linealpha=0.9,
        color = su4_color,
    )
    line_qu = (;
        label = "Quantum",
        linewidth=1.00,
        linealpha=1.0,
        color = :red,
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
    plot!([0., 1], [100., 100]; line_cl_label...) # Just so legend has no transparency
    plot!(kTs_ref, E_fm_cl; ylabel = L"$\left\langle E \right\rangle$", line_cl...)
    plot!(kTs_ref, E_fm_en; line_en...)
    plot!(kTs_ref, E_fm_qu; line_qu...)
    plot!(kTs, μs_fm_cl; yerr=sems_fm_cl, markerstrokecolor=dipole_color, error_bars...)
    plot!(kTs, μs_fm_en; yerr=sems_fm_en, markerstrokecolor=su4_color, error_bars...)

    # AFM (J > 0)
    p2 = plot(;
        yticks = (yticks_afm, [L"%$y" for y ∈ yticks_afm]),
        ylims = (-0.75*1.05, 0.05),
        xticks = (xticks, [L"%$x" for x ∈ xticks]),
        legend=:bottomright,
        params...
    )
    plot!([0., 1], [100., 100]; line_cl_label...) # Just so legend has no transparency
    plot!(kTs_ref, E_afm_cl; xlabel = L"k_bT", ylabel = L"$\left\langle E \right\rangle$", line_cl...)
    plot!(kTs_ref, E_afm_en; line_en...)
    plot!(kTs_ref, E_afm_qu; line_qu...)
    plot!(kTs, μs_afm_cl; yerr=sems_afm_cl, markerstrokecolor=dipole_color, error_bars...)
    plot!(kTs, μs_afm_en; yerr=sems_afm_en, markerstrokecolor=su4_color, error_bars...)

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
includet("SchrodingerLangevin.jl")
using .SchrodingerLangevin

using Random
using Plots
using LinearAlgebra
using FFTW
# using GLMakie 
using Interpolations
using Statistics

rng = MersenneTwister(111)

# Parameters
begin
    J = 1.0
    Δt = 0.01
    dur = 1000.0
    dur_bi = 50.0
    kTs = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
end

begin
    Es_all = []
    for kT ∈ kTs
        sys_e = entangled_pair(; J)
        rand!(sys_e)
        @time energy_trajectory!(sys_e, dur_bi, Δt, kT)
        @time Es = energy_trajectory!(sys_e, dur, Δt, kT)
        push!(Es_all, Es)
    end
end

begin
    ts = collect(0:length(Es_all[1])-1) .* Δt
    i_max = findfirst(x -> x > 200, ts)
    fig = GLMakie.Figure(resolution=(800,600))
    num_cols = 2
    for (i, Es) ∈ enumerate(Es_all)
        ax = GLMakie.Axis(fig[fldmod1(i, num_cols)...])
        es = fft(Es)
        ac = real(ifft(es .* conj.(es)))
        lines!(ax, ts[1:i_max], ac[1:i_max])
    end
    fig
end

begin
    J = 1.0
    kTs = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    # kTs = 0.1:0.1:1.0 
    bin_dur = [40.0, 25.0, 10.0, 6.0, 5.0, 5.0, 5.0]
    itp = LinearInterpolation(kTs, bin_dur)
    bin_func = kT -> kT < 0.1 ? 1.0 : itp(kT)
    dur_bi = 50.0
    Δt = 0.1
    num_samples = 1000
    kTs = 0.0:0.1:2.0
end

begin
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
end


begin
    kTs_ref = 0.000:0.001:2.0

    ## FM
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


    # scatter(kTs, μs)
    # errorbars!(kTs, μs, σs)
    layout = @layout [a; b]
    p = plot(p1, p2;
        layout,
        title = ["FM" "AFM"],
        figsize = (100,200)
    )
end

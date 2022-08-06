using GLMakie
using Sunny
import Statistics: mean
using LinearAlgebra
using JLD2
using Parameters

includet("model.jl")
includet("sim_utils.jl")
includet("plaquette_utils.jl")
includet("plot_utils.jl")



function bin_idcs(Zs, idx; width = 10)
    num_steps = size(Zs)[end]
    i1 = max(1, idx - round(Int, width/2) + 1)
    i2 = min(num_steps, i1 + width)
    i1:i2
end

# Generate and save trajectory
begin
    # dur = 500.0
    # dims = (100,100,1)
    dur = 100.0
    dims = (20,20,1)
    h = 15.35
    Δt = 0.01
    α = 0.1

    #= Set up system and initialize =#
    seed = 111
    rng = MersenneTwister(seed)
    sys = su3_skyrmion_model(dims; h, rng)
    rand!(sys)
end

begin
    #= Run and save trajectory =#
    integrator = LangevinHeunP(sys, kT, α)
    Zs = ket_trajectory!(integrator, Δt, dur)

    #= Save results =#
    save("trajectory_small.jld2", "Zs", Zs, "Δt", Δt, "h", h)
end

begin
    data = load("trajectory.jld2")
    @unpack Zs, Δt, h = data
    Sunny.init_from_coherents!(sys, Zs[:,:,:,:,end])
end


begin
    τ = 4.00
    base = 2.0
    # τs = [τ*(base^k) for k ∈ 0:2]   
    τs = [√10, 10.0, 100]
    idcs = [round(Int64, τ/Δt) for τ ∈ τs]

    Z_frames = [reshape(mean(Zs[:,:,:,:,bin_idcs(Zs, i; width=10)], dims=5), size(Zs)[1:4])
        # for i ∈ [1; idcs; 50000]
        for i ∈ idcs
    ]
    # push!(Z_frames, Zs[:,:,:,:,end])
    for frame ∈ Z_frames
        χ = plaquette_map(berry, frame)
        println("χ: ", sum(χ))
        println("Charge count: ", sum(χ)/2π)
    end
    vecs = sys.lattice.lat_vecs
    v₁, v₂ = vecs[:,1], vecs[:,2]
    fig, ax = plot_chirality_multi(Z_frames, v₁, v₂; offset=10, resolution=(1200,300), colorscheme=ColorSchemes.RdBu, clims=(-1.0, 1.0))

    # scene = LScene(fig[1,end+1]; show_axis=false)
    # plot_spins_color!(scene, Zs[:,:,:,:,end], sys; colorscheme=ColorSchemes.RdBu, arrowsize=0.5)
    fig
end

begin
    fig2 = Figure()
    ax2 = LScene(fig2[1,1]; show_axis=false)
    plot_spin_fluctuations!(ax2, Zs[:,:,:,:,end], sys;
        colorscheme=ColorSchemes.RdBu,
        tensor_scale=0.4,
        arrowlength=1.0,
        linewidth=0.08,
        arrowsize=0.3,
    )
    fig2
end
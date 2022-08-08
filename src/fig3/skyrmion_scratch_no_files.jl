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
    τ = 4.00
    base = 2.0
    τs = [4, 16, 256.0] 
    idcs = [round(Int64, τ/Δt) for τ ∈ τs]
    dims = (100,100,1)
    Δt = 0.01
    h = 15.35
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
    t = 0.0
    frames = []
    for i ∈ 1:length(τs)
        dur = i == 1 ? τs[1] : τs[i] - sum(τs[1:i-1])
        num_steps = round(Int, dur/Δt)
        for _ ∈ 1:num_steps
            evolve!(integrator, Δt)
        end
        push!(frames, sys._coherents)
    end

end

#= Save results =#
save("frames_4_16_256.jld2", "frames", frames, "Δt", Δt, "h", h)

begin
    data = load("frames_4_16_256.jld2")
    @unpack frames, Δt, h = data
    Sunny.init_from_coherents!(sys, frames[end])
end


begin
    τ = 4.00
    base = 2.0
    # τs = [τ*(base^k) for k ∈ 0:2]   
    τs = [4, 16, 256.0] 
    # τs = [3, 10, 100.0] 
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
    factor = 10 
    fig, ax = plot_chirality_multi(Z_frames, v₁, v₂; offset=10, resolution=(1200*factor,300*factor), colorscheme=ColorSchemes.RdBu, clims=(-1.0, 1.0))
    save("skyrmion_panels.png", fig)

    resolution = (round(Int, 400*factor), round(Int, 346.4101617*factor))
    fig1 = plot_chirality(Z_frames[1], v₁, v₂; resolution)
    save("skyrmion_panel_1.png", fig1)
    fig2 = plot_chirality(Z_frames[2], v₁, v₂; resolution)
    save("skyrmion_panel_2.png", fig2)
    fig3 = plot_chirality(Z_frames[3], v₁, v₂; resolution)
    save("skyrmion_panel_3.png", fig3)


    # scene = Axis(fig[1,end+1])
    # plot_spins_color!(scene, Zs[:,:,:,:,end], sys; colorscheme=ColorSchemes.RdBu, arrowsize=0.5)
    fig
end

begin
    Z = Zs[:,:,:,:,25600]
    fig = Figure()
    fig, ax = plot_chirality(Z, sys)
    ax = Axis(fig[2,1], aspect=1+sin(π/3))
    hidespines!(ax)
    hidedecorations!(ax)
    plot_spins_color!(ax, Z, sys; colorscheme=ColorSchemes.RdBu, arrowsize=0.5)
end

begin
    Z = Zs[:,:,:,:,25600]
    fig = Figure()
    fig, ax = plot_chirality(Z, sys; aspect=nothing, clims=(-1, 1))
    # ax = Axis(fig[2,1], aspect=1+sin(π/3))
    # hidespines!(ax)
    # hidedecorations!(ax)
    plot_spins_color!(ax, Z, sys; colorscheme=ColorSchemes.RdBu,
        arrowsize=0.4,
        linewidth=0.1,
        arrowlength=1.0,
    )

end

begin
    cam = cameracontrols(ax.scene)
    x, y, z = 115, 63, 40 # 4 blue
    # x, y, z = 121.5, 65, 20 # single blue bushel of wheat
    # x, y, z = 120, 81.37, 20 # single blue perfect spiral
    # x, y, z = 59, 44.00, 45 # half/half metastable structure
    # x, y, z = 66.75, 30.50, 30 # single opening 
    cam.lookat[] = [x, y, 0.0]
    cam.eyeposition[] = [x, y, z]
    update_cam!(ax.scene, cam)
end
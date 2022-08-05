# module 

using GLMakie
using Sunny
import Statistics: mean

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

# function fig3_quench(; dims=(100,100,1), kT=0.00, Δt=0.01, α=0.1, τ=0.3, base=5.0, num_frames=3, h=15.35, seed=111)
begin
    dims=(100,100,1)
    kT=0.00
    Δt=0.01
    α=0.1
    τ=0.1
    base=10.0
    num_frames=3
    h=15.35
    seed=111

    #= Frame calculations =#
    τs = [τ*base^k for k ∈ 1:num_frames]
    τ_idcs = [round(Int, τ/Δt) for τ ∈ τs]
    dur = τs[end]
    println("Duration: ", dur)

    #= Set up system and initialize =#
    rng = MersenneTwister(seed)
    sys = su3_skyrmion_model(dims; h, rng, rescaling=2.0)
    rand!(sys)

    #= Run and save trajectory =#
    integrator = LangevinHeunP(sys, kT, α)
    @time Zs = ket_trajectory!(integrator, Δt, dur)

    #= Choose frames from trajectory =#
    vecs = sys.lattice.lat_vecs
    v₁, v₂ = vecs[:,1], vecs[:,2]

    #= Plot Berry curvature per plaquette =#
    Z_frames = [reshape(mean(Zs[:,:,:,:,bin_idcs(Zs, i; width=10)], dims=5), size(Zs)[1:4])
        for i ∈ τ_idcs
    ]
    for frame ∈ Z_frames
        χ = plaquette_map(berry, frame)
        println("χ: ", sum(χ))
        println("Charge count: ", sum(χ)/2π)
    end
    # Z_frames = [Zs[:,:,:,:,i]
    #     for i ∈ τ_idcs
    # ]
    # push!(Z_frames, Zs[:,:,:,:,end])
    fig = plot_chirality_multi(Z_frames, v₁, v₂; offset=10, resolution=(1200,300), colorscheme=ColorSchemes.RdBu, clims=(-1.0, 1.0))
end


function fig3_steps(;
    dims=(100,100,1),
    kTs = [(0.5)^k for k ∈ 2:8],
    num_avgs = 10,
)
    #= Set up system and initialize =#
    rng = MersenneTwister(111)
    sys = su3_skyrmion_model(dims; h=15.35, rng)
    rand!(sys)


    #= Run and save trajectory =#
    frames = []
    @time for (i, kT) ∈ enumerate(kTs)
        println("Trajectory $i of $(length(kTs))")
        Δt = 0.01
        α = 0.1
        durs = [10.0, 100.0, 100.0]
        integrator = LangevinHeunP(sys, kT, α)
        Zs = ket_trajectory!(integrator, Δt, durs[i]) # don't need to save trajectory here
        push!(frames, Zs[:,:,:,:,end])
        χ = plaquette_map(berry, frames[end])
        println("χ: ", sum(χ))
        println("Charge count: ", sum(χ)/2π)
    end

    #= Choose frames from trajectory =#
    vecs = sys.lattice.lat_vecs
    v₁, v₂ = vecs[:,1], vecs[:,2]

    #= Plot Berry curvature per plaquette =#
    plot_chirality_multi(frames, v₁, v₂; offset=10, resolution=(1200,300), colorscheme=ColorSchemes.RdBu, clims=(-1.0, 1.0))
end





begin
    dims=(100,100,1)
    kT=0.00
    Δt=0.01
    α=0.1
    τ=0.1
    base=10.0
    # h=15.35
    rescaling =  2.0    # 2.0, 1.6
    h = 15.35 
    D = 19.0 
    seed=111

    #= Frame calculations =#
    dur = 40.0

    #= Set up system and initialize =#
    rng = MersenneTwister(seed)
    sys = su3_skyrmion_model(dims; h, rng, rescaling)
    rand!(sys)

    #= Run and save trajectory =#
    integrator = LangevinHeunP(sys, kT, α)
    @time Zs = ket_trajectory!(integrator, Δt, dur)

    #= Choose frames from trajectory =#
    vecs = sys.lattice.lat_vecs
    v₁, v₂ = vecs[:,1], vecs[:,2]

    χ = plaquette_map(berry, Zs[:,:,:,:,end])
    println("χ: ", sum(χ))
    println("Charge count: ", sum(χ)/2π)

    #= Plot Berry curvature per plaquette =#
    fig = plot_chirality(Zs[:,:,:,:,end], v₁, v₂; colorscheme=ColorSchemes.RdBu, clims=(-1.0, 1.0))
end
module 

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

function fig2_quench(; dims=(100,100,1), kT=0.01, dur=100.0, Δt=0.01, α=0.1)
    #= Set up system and initialize =#
    rng = MersenneTwister(111)
    sys = su3_skyrmion_model(dims; h=15.25, rng)
    rand!(sys)

    #= Run and save trajectory =#
    integrator = LangevinHeunP(sys, kT, α)
    Zs = ket_trajectory!(integrator, Δt, dur)

    #= Choose frames from trajectory =#
    vecs = sys.lattice.lat_vecs
    v₁, v₂ = vecs[:,1], vecs[:,2]
    τ = 0.3 
    τs = [τ*4^k for k ∈ 1:4]
    τ_idcs = [round(Int, τ/Δt) for τ ∈ τs]
    println(τ_idcs)
    println(size(Zs))

    #= Plot Berry curvature per plaquette =#
    Z_frames = [reshape(mean(Zs[:,:,:,:,bin_idcs(Zs, i; width=10)], dims=5), size(Zs)[1:4])
        for i ∈ τ_idcs
    ]
    # Z_frames = [Zs[:,:,:,:,i]
    #     for i ∈ τ_idcs
    # ]
    push!(Z_frames, Zs[:,:,:,:,end])
    plot_chirality_multi(Z_frames, v₁, v₂; offset=10, resolution=(1300,300))
end


function fig2_steps(;
    dims=(100,100,1),
    kTs = [(0.5)^k for k ∈ 2:8],
)
    #= Set up system and initialize =#
    rng = MersenneTwister(111)
    sys = su3_skyrmion_model(dims; h=15.25, rng)
    rand!(sys)


    #= Run and save trajectory =#
    frames = []
    @time for (i, kT) ∈ enumerate(kTs)
        println("Trajectory $i of $(length(kTs))")
        Δt = 0.01
        α = 0.1
        dur = 20.0 
        integrator = LangevinHeunP(sys, kT, α)
        Zs = ket_trajectory!(integrator, Δt, dur) # don't need to save trajectory here
        push!(frames, Zs[:,:,:,:,end])
    end

    #= Choose frames from trajectory =#
    vecs = sys.lattice.lat_vecs
    v₁, v₂ = vecs[:,1], vecs[:,2]

    #= Plot Berry curvature per plaquette =#
    plot_chirality_multi(frames, v₁, v₂; offset=10, resolution=(1300,300))
end

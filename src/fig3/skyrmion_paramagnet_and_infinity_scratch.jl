using GLMakie
using Sunny
import Statistics: mean
using LinearAlgebra
using JLD2
using Parameters
using Formatting
using DelimitedFiles

includet("model.jl")
includet("sim_utils.jl")
includet("plaquette_utils.jl")
includet("plot_utils.jl")



raw_angles = readdlm("src/fig3/hao_angles.dat")

function angle_to_cs(angles)
    θ, ϕ, α1, α2 = angles
    return SA[
        exp(im*α1)*sin(θ)*cos(ϕ),
        cos(θ),
        exp(im*α2)*sin(θ)*sin(ϕ)
]
end

function reinterpret_hao(raw_angles)
    angles = reinterpret(SVector{4, Float64}, raw_angles)
    angles_new = zeros(SVector{4, Float64}, 10, 10, 1, 1)
    for (i, angle) ∈ enumerate(angles)
        r = fld1(i, 10)
        c = mod1(i, 10)
        angles_new[c, r, 1, 1] = angles[i]
    end
    cs = map(angle_to_cs, angles_new)
end

function periodically_extend(gs, factor) 
    x, y = size(gs)[1:2]
    Zs = zeros(Sunny.CVec{3}, x*factor, y*factor, 1, 1)
    for i ∈ 1:factor, j ∈ 1:factor
        Zs[(i-1)*x+1:i*x, (j-1)*y+1:j*y, 1, 1] .= gs[:,:,1,1]
    end
    return Zs
end 

begin
    Zs = reinterpret_hao(raw_angles)
    Zs = periodically_extend(Zs, 10)
    dims_hao = (100, 100, 1)
    sys_hao = su3_skyrmion_model(dims_hao; D=19.0, h=15.349710144927537, rng)
    Sunny.init_from_coherents!(sys_hao, Zs)
    energy(sys_hao)/10000
end

function bin_idcs(Zs, idx; width = 10)
    num_steps = size(Zs)[end]
    i1 = max(1, idx - round(Int, width/2) + 1)
    i2 = min(num_steps, i1 + width)
    i1:i2
end

# Generate and save trajectory

begin
    trial = "hao"
    # trial = "qpm"
    # trial = "inf"

    dur = 500.0
    kT = 0.11
    dims = (100,100,1)
    h = 15.35
    # h = 15.20
    Δt = 0.01
    α = 0.1

    #= Set up system and initialize =#
    seed = 111
    rng = MersenneTwister(seed)
    sys = su3_skyrmion_model(dims; h, rng)
    if trial == "qpm"
        Z₀ = fill(SVector{3, ComplexF64}([0, 1, 0]), size(sys._coherents))
        Sunny.init_from_coherents!(sys, Z₀)
    elseif trial == "hao"
        raw_angles = readdlm("src/fig3/hao_angles.dat")
        Z₀ = reinterpret_hao(raw_angles)
        Z₀ = periodically_extend(Z₀, 10)
        println(size(Z₀))
        Sunny.init_from_coherents!(sys, Z₀)
    elseif trial == "inf"
        rand!(sys)
    else
        println("What're you doing, man?")
    end
end

begin
    num_steps = round(Int, dur/Δt)
    quench_steps = round(Int, 10.0/Δt)
    integrator = LangevinHeunP(sys, kT, α)

    trajectory = zeros(SVector{3, ComplexF64}, size(sys._coherents)..., num_steps + quench_steps)
    kTs_all = zeros(num_steps + quench_steps)

    integrator.kT = kT
    @time for c ∈ 1:(num_steps + quench_steps)
        kT′ = c < num_steps + 1 ? kT : 0.0
        integrator.kT = kT′
        evolve!(integrator, Δt)
        trajectory[:,:,:,:,c] .= sys._coherents
        kTs_all[c] = kT′
    end
end

begin
    χs = zeros(size(trajectory)[end])
    for i ∈ 1:size(trajectory)[end]
        χ = plaquette_map(berry, trajectory[:,:,:,:,i])
        χs[i] = round(Int, sum(χ)/2π)
    end
end


begin
    text = ["Total charge: $(χ)\nkT=$(format("{:.4f}", T))" for (χ, T) ∈ zip(χs, kTs_all)]
    lat_vecs = sys.lattice.lat_vecs
    v1, v2 = lat_vecs[:,1], lat_vecs[:,2]
    animate_chirality(trajectory, v1, v2;
        filename="$(trial)_kT=$kT.mp4", clims=(-1, 1), skip_interval=100,
        text
    )
end

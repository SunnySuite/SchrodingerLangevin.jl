using Sunny
using GLMakie
using LinearAlgebra
using JLD2 
using Formatting
using FFTW


includet("model.jl")
includet("sim_utils.jl")
includet("plaquette_utils.jl")
includet("plot_utils.jl")

################################################################################
# Relax into ground state
################################################################################

begin
    kTs = [10.0 * 0.85 ^ k for k ∈ 0:50]
    push!(kTs, 0.0)
    step_dur = 10.0
    dims = (100,100,1)
    h = 15.35
    D = 19.0
    Δt = 0.01
    α = 0.1

    #= Set up system and initialize =#
    seed = 111
    rng = MersenneTwister(seed)
    sys = su3_skyrmion_model(dims; D, h, rng)
    rand!(sys)
    kTs
end

begin
    num_steps = round(Int, step_dur/Δt)
    integrator = LangevinHeunP(sys, kT, α)

    for kT ∈ kTs
        println("kT = $kT...")
        integrator.kT = kT
        @time for c ∈ 1:num_steps
            evolve!(integrator, Δt)
        end
    end
end

begin
    num_steps = round(Int, 100.0/Δt)
    integrator.kT = 0.0
    @time for c ∈ 1:num_steps
        evolve!(integrator, Δt)
    end

    mags = norm.(sys._dipoles)
    μ_dipole =  mean(mags)
end

plot_spins(sys; arrowlength=1.0)

gs_langevin = copy(sys._coherents)
energy(sys)/10000


################################################################################
# Hao GS
################################################################################
using DelimitedFiles

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

begin
    Zs = reinterpret_hao(raw_angles)
    dims_hao = (10, 10, 1)
    sys_hao = su3_skyrmion_model(dims_hao; D=19.0, h=15.349710144927537, rng)
    Sunny.init_from_coherents!(sys_hao, Zs)

    M = sum(sys_hao._dipoles)
    M = norm(M)/100
    E = energy(sys_hao)/100.0
    println("M=$M, E=$E")
end

plot_spins(sys_hao; arrowlength=0.5)

# Relax
begin
    integrator = LangevinHeunP(sys_hao, 0.0, 0.1)
    dur = 100.0
    num_steps = round(Int, dur/Δt)
    @time for _ ∈ 1:num_steps
        evolve!(integrator, Δt)
    end
end
energy(sys_hao)/100
plot_spins(sys_hao; arrowlength=0.5)



################################################################################
# Scale up hao
################################################################################
raw_angles = readdlm("src/fig3/hao_angles.dat")


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

begin
    integrator = LangevinHeunP(sys_hao, 0.0, 0.1)
    dur = 100.0
    num_steps = round(Int, dur/Δt)
    @time for _ ∈ 1:num_steps
        evolve!(integrator, Δt)
    end
    energy(sys_hao)/10000
end

plot_spins(sys_hao; arrowlength=0.5)

begin
    Zs = sys_hao._coherents
    vecs = sys_hao.lattice.lat_vecs
    v1, v2 = vecs[:,1], vecs[:,2]
    fig, ax = plot_chirality(Zs, v1, v2; clims=(-0.3, 0.3))

    # plot_spins_color!(ax, Zs, sys_hao)
    fig
end

dipoles = sys_hao._dipoles
xs = [s[1] for s in dipoles]
ys = [s[2] for s in dipoles]
zs = [s[3] for s in dipoles]

Xs = fft(xs)
Ys = fft(ys)
Zs = fft(zs)

S = Xs .* conj.(Xs) + Ys .* conj.(Ys) + Zs .* conj.(Zs)

heatmap(S[:,:,1,1])
using DrWatson
@quickactivate "SU3SkyrmionFormation"

using Sunny


function anneal!(sampler, kTs)
    Es = zeros(length(kTs))        
    for (i, kT) in enumerate(kTs)
          set_temp!(sampler, kT)
          sample!(sampler)                    
          Es[i] = running_energy(sampler)   
    end
    return Es    
end

function dipole_trajectory!(integrator, Δt, dur)
    sys = integrator.sys
    n = round(Int, dur/Δt)
    trajectory = zeros(Sunny.Vec3, size(sys)..., n)
    for i ∈ 1:n
        evolve!(integrator, Δt)
        trajectory[:,:,:,:,i] .= sys._dipoles 
    end
    trajectory
end

function ket_trajectory!(integrator, Δt, dur)
    sys = integrator.sys
    n = round(Int, dur/Δt)
    trajectory = zeros(eltype(sys._coherents), size(sys)..., n)
    for i ∈ 1:n
        evolve!(integrator, Δt)
        trajectory[:,:,:,:,i] .= sys._coherents
    end
    trajectory
end

function energy_trajectory!(integrator, Δt, dur)
    sys = integrator.sys
    n = round(Int, dur/Δt)
    Es = zeros(n)
    for i ∈ 1:n
        evolve!(integrator, Δt)
        Es[i] = energy(sys)
    end
    Es
end

function thermalize!(integrator, Δt, dur, kT)
    integrator.kT = kT  # write set_temp! function?
    n = round(Int, dur/Δt)
    for _ ∈ 1:n
        evolve!(integrator, Δt)
    end
    nothing
end
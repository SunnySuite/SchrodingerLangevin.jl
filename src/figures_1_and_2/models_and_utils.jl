################################################################################
# Functions for generating statistics
################################################################################
function energy_trajectory!(sys, dur, Δt, kT)
    n = round(Int, dur/Δt)
    Es = zeros(n+1)
    Es[1] = energy(sys)
    for i ∈ 1:n
        heun_langevin_step!(sys, Δt, kT)
        Es[i+1] = energy(sys)
    end
    return Es
end

function sample_variance(vals)
    N = length(vals)
    μ = sum(vals) / N 
    (1/(N-1)) * sum((vals .- μ) .^ 2)
end

function generate_statistics(Δt, num_samples, kTs;
    sys_func, dur_trajectory = 10.0, dur_burnin=10.0,
)
    μs = zero(kTs)
    sems = zero(kTs)
    for (i, kT) in enumerate(kTs)
        println("Collecting statistics for kT=$kT...")

        sys = sys_func()
        rand!(sys)
        energy_trajectory!(sys, dur_burnin, Δt, kT)

        Es = zeros(num_samples)
        for j ∈ 1:num_samples
            E = energy_trajectory!(sys, dur_trajectory, Δt, kT)
            Es[j] = mean(E) 
        end
        μs[i] = mean(Es)
        sems[i] = √(sample_variance(Es)/num_samples)
    end
    return (; μs, sems)
end

################################################################################
# Models
################################################################################
function classical_pair(; J=1.0, α=0.1, rng=nothing)
    System(; N=2, L=2, J, rng, α)
end

function entangled_pair(; J=1.0, α=0.1, rng=nothing)
    v = [0.0 1.0 -1.0 0.0] / √2
    Λ = (J/4)*I(4) - J*v'*v
    System(; N=4, L=1, J=0.0, Λ=[Λ], rng, α)
end

function dipole_anisotropy(; D=-1.0, α=0.1, rng=nothing)
    System(; N=2, L=1, D, rng, spin_rescaling=2.0, α)
end

function su3_anisotropy(; D=-1.0, α=0.1, rng=nothing)
    Λ = zeros(ComplexF64, 3, 3)
    Λ[1,1] = Λ[3,3] = D  # D*(Sᶻ)² in spin-1 representation
    System(; N=3, L=1, Λ=[Λ], D=0.0, rng, α)
end


################################################################################
# Analytical energy functions
################################################################################
function energy_pair_dipole(β, J=1.0)
    (1/β) - (J/4)*coth(β*J/4)
end

function energy_pair_entangled(β, J=1.0)
    if β >= 1000 || isinf(β)
        return J > 0 ? -0.75 : -0.25
    end
    num = 24 + 18J*β + 6(J*β)^2 + (J*β)^3 + 6*exp(J*β)*(J*β - 4)
    denom = 4*β*(2 - 2*exp(J*β) + 2*J*β + (J*β)^2) 
    return num/denom
end

function energy_pair_quantum(β, J=1.0)
    if β >= 1000 || isinf(β)
        return J > 0 ? -0.75 : -0.25
    end
    -3J*(exp(β*J) - 1) / (4*(exp(β*J) + 3))
end

function energy_aniso_dipole(β, D=-1.0)
    if β >= 700 || isinf(β)
        return -1.0
    end
    real(1/(2β) - (√complex(D) * exp(-D*β))/(√(π*β)*erf(√(complex(D)*β))))
end

function energy_aniso_su3(β, D=-1.0)
    isinf(β) && (return -1.0)
    2/β + (D^2 * β)/(1 - exp(D*β) + D*β)
end

function energy_aniso_quantum(β, D=-1.0)
    isinf(β) && (return -1.0)
    2D/(2 + exp(D*β))
end
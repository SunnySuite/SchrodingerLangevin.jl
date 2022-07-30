include("SchrodingerLangevin.jl")
using .SchrodingerLangevin
using LinearAlgebra

function entangled_exchange(J)
    v = [0.0 1.0 -1.0 0.0] / √2
    return (J/4)*I(4) - J*v'*v
end

J = 1.0
Λ = entangled_exchange(J)

sys_c = System(; N=2, L=2, J)
sys_e = System(; N=4, L=1, J=0.0, Λ=[Λ])

Δt = 0.01
kT = 0.0

begin
    SchrodingerLangevin.heun_langevin_step!(sys_e, Δt, kT)
    SchrodingerLangevin.energy(sys_e)
end
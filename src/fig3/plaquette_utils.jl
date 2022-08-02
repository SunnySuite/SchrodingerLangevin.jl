function arg(z::ComplexF64)
    x = real(z)
    y = imag(z)
    argument = 0.0
    if x > 0.0
        argument = atan(y/x)
    elseif x <  0.0 && y >= 0
        argument = atan(y/x) + π
    elseif x < 0.0 && y < 0
        argument = atan(y/x) - π
    elseif x == 0.0 && y > 0
        argument = π/2
    elseif x == 0.0 && y < 0
        argument = -π/2
    else
        argument = NaN
    end
    argument
end


"""
    berry(o₁, o₂, o₃)
Determine the Berry curvature for three oriented SU(3) spins.
"""
function berry(o₁, o₂, o₃)
   n₁ = o₁' * o₂
   n₂ = o₂' * o₃
   n₃ = o₃' * o₁
   2 * arg(n₁ * n₂ * n₃) 
end


function plaquette_idcs(dims::Tuple{Int, Int, Int})
    dx, dy, dz = dims
    (dz != 1) && println("Warning: Ignoring lattice c vector.")
    Triple = Tuple{Int, Int, Int}
    idcs = Array{Tuple{Triple, Triple, Triple}, 4}(undef, 2, dx+1, dy+1, 1) 
    c = 1 
    for j ∈ 1:(dy+1)
        for i ∈ 1:(dx+1)
            idcs[1, i, j, 1] = (
                (mod1(i, dx),   mod1(j, dy),   1),
                (mod1(i+1, dx), mod1(j, dy),   1),
                (mod1(i, dx),   mod1(j+1, dy), 1),
            )
            idcs[2, i, j, 1] = (
                (mod1(i+1, dx), mod1(j, dy),   1),
                (mod1(i+1, dx), mod1(j+1, dy), 1),
                (mod1(i,   dx), mod1(j+1, dy), 1),
            )
        end
    end
    return idcs
end

plaquette_idcs(x::Array{T, 4}) where T= plaquette_idcs(size(x)[1:3]) 


function plaquette_map(f::Function, x::Array{T, 4}) where T
    dims = size(x)
    @assert dims[3] == dims[4] == 1 "Multiple sites and multiple layers are not supported."
    dims = dims[1:3]
    out = zeros(Float64, 2, dims...) # How to set type based on return type of f?
    idcs_all = plaquette_idcs(x)
    for i ∈ CartesianIndices(dims) 
        idcs = idcs_all[1, i]
        out[1,i] = f(x[idcs[1]...,1], x[idcs[2]...,1], x[idcs[3]...,1])
        idcs = idcs_all[2, i]
        out[2,i] = f(x[idcs[1]...,1], x[idcs[2]...,1], x[idcs[3]...,1])
    end
    out
end
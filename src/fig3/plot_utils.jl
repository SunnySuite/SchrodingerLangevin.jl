using GLMakie
using Observables
using ColorSchemes
using ColorTypes
import Makie: get_dim, surface_normals
using GeometryBasics
using StaticArrays

function plot_chirality(Z, v₁, v₂;
    colorscheme=ColorSchemes.RdBu,
    resolution=(400, 400),
    aspect = nothing,
    clims = nothing,
)
    dims = size(Z)
    nx, ny = dims[1:2]
    v₁ = Point3f(v₁)
    v₂ = Point3f(v₂)
    plaq1(p) = GLMakie.Polygon(Point2f.([p, p+v₁, p+v₂]))
    plaq2(p) = GLMakie.Polygon(Point2f.([p+v₁, p+v₁+v₂, p+v₂]))
    
    χ = plaquette_map(berry, Z)

    if isnothing(clims)
        max = maximum(abs.(χ))
        clims = (-max, max)
    end

    pgons = GLMakie.Polygon[]
    colors = ColorTypes.RGB{Float64}[]
    for r ∈ 1:nx
        for c ∈ 1:ny
            base = (r-1)*v₁ + (c-1)*v₂
            push!(pgons, plaq1(base))
            push!(colors, get(colorscheme, χ[1,r,c,1,1], clims))
            push!(pgons, plaq2(base))
            push!(colors, get(colorscheme, χ[2,r,c,1,1], clims))
        end
    end

    base = (0, 0) 
    corner = (nx-1)*v₁ + (ny-1)*v₂
    x1, x2 = base[1], corner[1]
    y1, y2 = base[2], corner[2]
    isnothing(aspect) && (aspect = (x2-x1)/(y2-y1))

    fig = Figure(; resolution)
    # ax = Axis(fig[1,1]; aspect)
    ax = LScene(fig[1,1]; show_axis=false)
    # xlims!(ax, x1, x2)
    # ylims!(ax, y1, y2)
    # hidespines!(ax); hidedecorations!(ax)
    poly!(ax, pgons; color=colors)

    fig, ax
end

function plot_chirality(sys::Sunny.SpinSystem; kwargs...)
    Z = sys._coherents
    lat_vecs = sys.lattice.lat_vecs
    plot_chirality(Z, lat_vecs[:,1], lat_vecs[:,2]; kwargs...)
end

function plot_chirality(Zs, sys::Sunny.SpinSystem; kwargs...)
    lat_vecs = sys.lattice.lat_vecs
    plot_chirality(Zs, lat_vecs[:,1], lat_vecs[:,2]; kwargs...)
end



function plot_chirality_multi(Zs, v₁, v₂;
    colorscheme=ColorSchemes.RdBu,
    clims = (-2π, 2π),
    offset = 1,
    kwargs...
)
    num_panels = length(Zs)
    Z = Zs[1]
    dims = size(Z)
    nx, ny = dims[1:2]
    v₁ = Point3f(v₁)
    v₂ = Point3f(v₂)
    v_offset = v₁ * (nx+offset)

    plaq1(p) = GLMakie.Polygon(Point2f.([p, p+v₁, p+v₂]))
    plaq2(p) = GLMakie.Polygon(Point2f.([p+v₁, p+v₁+v₂, p+v₂]))

    base = (0, 0) 
    corner = (nx-1)*num_panels*v₁ + (num_panels-1)*offset*v₁ + (ny-1)*v₂
    x1, x2 = base[1], corner[1]
    y1, y2 = base[2], corner[2]
    aspect = (x2-x1)/(y2-y1)

    fig = Figure(; kwargs...)
    ax = Axis(fig[1,1:length(Zs)]; aspect)
    hidespines!(ax); hidedecorations!(ax)

    for (i, Z) ∈ enumerate(Zs)
        v₀ = (i-1) * v_offset
        Χ = (plaquette_map(berry, Z))
        # max = maximum(abs.(Χ))
        # clims = (-max, max)
        pgons = GLMakie.Polygon[]
        colors = ColorTypes.RGB{Float64}[]
        for r ∈ 1:nx
            for c ∈ 1:ny
                base = (r-1)*v₁ + (c-1)*v₂ + v₀
                push!(pgons, plaq1(base))
                push!(colors, get(colorscheme, Χ[1,r,c,1,1], clims))
                push!(pgons, plaq2(base))
                push!(colors, get(colorscheme, Χ[2,r,c,1,1], clims))
            end
        end
        poly!(ax, pgons; color=colors)
    end
    # Colorbar(fig[1,2]; colormap=colorscheme, colorrange=clims)
    fig, ax
end


function get_colors(colorscheme, x)
    nx, ny = size(x)[2:3]
    colors = Array{ColorTypes.RGB{Float64}, 1}(undef, nx*ny*2)
    count = 1
    for c ∈ 1:ny, r ∈ 1:nx
        colors[count] = get(colorscheme, x[1,r,c,1,1])
        count += 1
        colors[count] = get(colorscheme, x[2,r,c,1,1])
        count += 1
    end
    colors
end

function animate_chirality(Zs, v₁, v₂;
    filename = "chirality_anim.mp4",
    colorscheme=ColorSchemes.viridis,
    framerate=30,
    skip_interval=1,
)
    dims = size(Zs)[1:4]
    nx, ny = dims[1:2]
    v₁ = Point3f(v₁)
    v₂ = Point3f(v₂)
    plaq1(p) = GLMakie.Polygon(Point2f.([p, p+v₁, p+v₂]))
    plaq2(p) = GLMakie.Polygon(Point2f.([p+v₁, p+v₁+v₂, p+v₂]))
    idx = Observable(1)
    Χ = @lift(plaquette_map(berry, Zs[:,:,:,:,$idx]))
    colors = @lift(get_colors(colorscheme, $Χ))
    pgons = GLMakie.Polygon[]
    for c ∈ 1:ny, r ∈ 1:nx
        base = (r-1)*v₁ + (c-1)*v₂
        push!(pgons, plaq1(base))
        push!(pgons, plaq2(base))
    end

    fig = Figure()
    ax = Axis(fig[1,1])
    hidespines!(ax); hidedecorations!(ax)
    poly!(ax, pgons; color=colors)

    numsteps = size(Zs)[end]
    GLMakie.record(fig, filename, 1:skip_interval:numsteps; framerate) do i
        idx[] = i
    end
end

function animate_chirality(Zs, sys; kwargs...)
    lat_vecs = sys.lattice.lat_vecs
    animate_chirality(Zs, lat_vecs[:,1], lat_vecs[:,2]; kwargs...)
end



begin
    dims = (100,100,1)
    rng = MersenneTwister(111)
    sys = su3_skyrmion_model(dims; h=15.35, rng)
    rand!(sys)
end

#= Run and save trajectory =#
begin
    dur = 40.0
    Δt = 0.004
    kT = 0.0
    integrator = LangevinHeunP(sys, kT, 0.1)
    Zs = ket_trajectory!(integrator, Δt, dur)
    Z = Zs[:,:,:,:,end]
end


function plot_spins_color(Zs, sys;
    resolution=(600,400),
    kwargs...
)
    fig = GLMakie.Figure(; resolution)
    ax = GLMakie.LScene(fig[1,1]; show_axis=false, kwargs...)

    plot_spins_color!(ax, Zs, sys; resolution, kwargs...)
    fig
end

function plot_spins_color!(ax, Zs, sys;
    colorscheme = ColorSchemes.diverging_linear_bjr_30_55_c53_n256,
    arrowlength=1.5,
    arrowsize=0.3,
    linewidth=0.15,
    kwargs...
)
    points = GLMakie.Point3f0.(vec(sys.lattice))
    vecs = GLMakie.Vec3f0.(vec(sys._dipoles))
    points .-= points[1]
    lengths = norm.(vecs)
    lengths .-= minimum(lengths)
    lengths ./= maximum(lengths)
    # color = get(ColorSchemes.viridis, lengths)
    # color = get(colorscheme, lengths)
    # color = [(get(colorscheme, 2*length), length) for length ∈ lengths]
    color = [(:firebrick2, length) for length ∈ lengths]
    GLMakie.arrows!(
        ax, points, vecs;
        color,
        linewidth, arrowsize,
        lengthscale=arrowlength,
        kwargs...,
    )
    fig
end


expectation(op, ψ) = real(ψ' * op * ψ)
comm(a,b)     = a*b-b*a
anticomm(a,b) = a*b + b*a

function create_mesh(x, y, z)
    positions = vec(map(CartesianIndices(z)) do i
    GeometryBasics.Point{3, Float32}(
        get_dim(x, i, 1, size(z)),
        get_dim(y, i, 2, size(z)),
        z[i])
    end)
    faces = decompose(GLTriangleFace, Rect2D(0f0, 0f0, 1f0, 1f0), size(z))
    normals = surface_normals(x, y, z)
    vertices = GeometryBasics.meta(positions; normals=normals)
    meshObj = GeometryBasics.Mesh(vertices, faces)
    meshObj
end

function T(ψ, sx, sy, sz)
    t11 = expectation(sx*sx, ψ) - expectation(sx, ψ)^2
    t12 = 0.5*expectation(anticomm(sx,sy), ψ) - expectation(sx, ψ)*expectation(sy, ψ)
    t13 = 0.5*expectation(anticomm(sx,sz), ψ) - expectation(sx, ψ)*expectation(sz, ψ)
    t22 = expectation(sy*sy, ψ) - expectation(sy, ψ)^2
    t23 = 0.5*expectation(anticomm(sz,sy), ψ) - expectation(sz, ψ)*expectation(sy, ψ)
    t33 = expectation(sz*sz, ψ) - expectation(sz, ψ)^2
    SMatrix{3,3,Float64,9}([t11 t12 t13;
                            t12 t22 t23;
                            t13 t23 t33])
end

function ellipsoid_point(t, θ, ϕ)
    t[:,1]*cos(θ)*cos(ϕ) + t[:,2]*cos(θ)*sin(ϕ) + t[:,3]*sin(θ) 
end

function points_from_tensor2(t, f0=SA[0.0, 0.0, 0.0]; num_θ=25, num_ϕ=50, scale=1.0)
    θs = LinRange(-π/2, π/2, num_θ)
    ϕs = LinRange(0.0, 2π, num_ϕ)
    f1, f2, f3 = t[:,1], t[:,2], t[:,3]

    points = [ellipsoid_point(t, θ, ϕ) for θ ∈ θs, ϕ ∈ ϕs] .* scale
    xs = [f0[1] + points[i,j][1] for i ∈ 1:num_θ, j ∈ 1:num_ϕ]
    ys = [f0[2] + points[i,j][2] for i ∈ 1:num_θ, j ∈ 1:num_ϕ]
    zs = [f0[3] + points[i,j][3] for i ∈ 1:num_θ, j ∈ 1:num_ϕ]
    create_mesh(xs, ys, zs)
end


function plot_spin_fluctuations!(ax, Zs, sys;
    colorscheme = ColorSchemes.diverging_linear_bjr_30_55_c53_n256,
    arrowlength=1.0,
    arrowsize=0.3,
    linewidth=0.15,
    tensor_scale=1.0,
    kwargs...
)
    points = GLMakie.Point3f0.(vec(sys.lattice))
    vecs = GLMakie.Vec3f0.(vec(sys._dipoles))
    lengths = norm.(vecs)
    lengths .-= minimum(lengths)
    lengths ./= maximum(lengths)
    # color = get(ColorSchemes.viridis, lengths)
    color = [(get(colorscheme, length), length)  for length ∈ lengths]
    GLMakie.arrows!(
        ax, points, vecs;
        color,
        linewidth, arrowsize,
        lengthscale=arrowlength,
        kwargs...,
    )

    S = Sunny.gen_spin_ops(3)
    dipoles = sys._dipoles
    for i ∈ eachindex(Zs)
        Z = Zs[i]
        t = T(Z, S[1], S[2], S[3])
        mesh_points = points_from_tensor2(t, arrowlength*dipoles[i] + points[i]; scale = tensor_scale)
        mesh!(ax, mesh_points, color=(color[i][1], color[i][2]*0.5), shading=true, transparency=true) 
    end

end
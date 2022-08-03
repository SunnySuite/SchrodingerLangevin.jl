using GLMakie
using Observables
using ColorSchemes
using ColorTypes

function plot_chirality(Z, v₁, v₂; colorscheme=ColorSchemes.viridis)
    dims = size(Z)
    nx, ny = dims[1:2]
    v₁ = Point3f(v₁)
    v₂ = Point3f(v₂)
    plaq1(p) = GLMakie.Polygon(Point2f.([p, p+v₁, p+v₂]))
    plaq2(p) = GLMakie.Polygon(Point2f.([p+v₁, p+v₁+v₂, p+v₂]))
    #idx = Node(1)
    Χ = plaquette_map(berry, Z)
    pgons = GLMakie.Polygon[]
    colors = ColorTypes.RGB{Float64}[]
    for r ∈ 1:nx
        for c ∈ 1:ny
            base = (r-1)*v₁ + (c-1)*v₂
            push!(pgons, plaq1(base))
            push!(colors, get(colorscheme, Χ[1,r,c,1,1]))
            push!(pgons, plaq2(base))
            push!(colors, get(colorscheme, Χ[2,r,c,1,1]))
        end
    end

    fig = Figure()
    ax = Axis(fig[1,1])
    hidespines!(ax); hidedecorations!(ax)
    poly!(ax, pgons; color=colors)

    fig
end


function plot_chirality_multi(Zs, v₁, v₂;
    colorscheme=ColorSchemes.viridis,
    offset = 1,
    kwargs...
)
    Z = Zs[1]
    dims = size(Z)
    nx, ny = dims[1:2]
    v₁ = Point3f(v₁)
    v₂ = Point3f(v₂)
    v_offset = v₁ * (nx+offset)

    plaq1(p) = GLMakie.Polygon(Point2f.([p, p+v₁, p+v₂]))
    plaq2(p) = GLMakie.Polygon(Point2f.([p+v₁, p+v₁+v₂, p+v₂]))

    fig = Figure(; kwargs...)
    ax = Axis(fig[1,1])
    hidespines!(ax); hidedecorations!(ax)

    for (i, Z) ∈ enumerate(Zs)
        v₀ = (i-1) * v_offset
        Χ = plaquette_map(berry, Z)
        pgons = GLMakie.Polygon[]
        colors = ColorTypes.RGB{Float64}[]
        for r ∈ 1:nx
            for c ∈ 1:ny
                base = (r-1)*v₁ + (c-1)*v₂ + v₀
                push!(pgons, plaq1(base))
                push!(colors, get(colorscheme, Χ[1,r,c,1,1]))
                push!(pgons, plaq2(base))
                push!(colors, get(colorscheme, Χ[2,r,c,1,1]))
            end
        end
        poly!(ax, pgons; color=colors)
    end

    fig
end


function plot_chirality(sys::Sunny.SpinSystem)
    Z = sys._coherents
    lat_vecs = sys.lattice.lat_vecs
    plot_chirality(Z, lat_vecs[:,1], lat_vecs[:,2])
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



dims = (20,20,1)
rng = MersenneTwister(111)
sys = su3_skyrmion_model(dims; h=15.25, rng)
rand!(sys)

#= Run and save trajectory =#
begin
    dur = 100.0
    Δt = 0.004
    kT = 0.0
    integrator = LangevinHeunP(sys, kT, 0.1)
    Zs = ket_trajectory!(integrator, Δt, dur)
    Z = Zs[:,:,:,:,end]
end


begin 
    function plot_spins_color(Zs, sys;
        resolution=(600,400),
        kwargs...
    )
        fig = GLMakie.Figure(; resolution)
        ax = GLMakie.LScene(fig[1,1]; show_axis=false, kwargs...)

        points = GLMakie.Point3f0.(vec(sys.lattice))
        vecs = GLMakie.Vec3f0.(vec(sys._dipoles))
        lengths = norm.(vecs)
        lengths .-= minimum(lengths)
        lengths ./= maximum(lengths)
        color = get(ColorSchemes.viridis, lengths)

        GLMakie.arrows!(
            ax, points, vecs;
            color,
        )
        fig
    end
    plot_spins_color(Z, sys)
end
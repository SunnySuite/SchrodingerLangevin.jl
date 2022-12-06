################################################################################
# Functions for calculating over plaquettes, in particular Berry curvature
################################################################################
function berry(o₁, o₂, o₃)
   n₁ = o₁' * o₂
   n₂ = o₂' * o₃
   n₃ = o₃' * o₁
   angle(n₁ * n₂ * n₃) 
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


################################################################################
# Plotting functions
################################################################################
function aspect_ratio(x_panel, y_panel, x_offset, y_offset, numrows, numcols)
    corners = [
        (0, 0) ,
        numrows*y_panel + (numrows-1)*y_offset,
        numcols*x_panel + (numcols-1)*x_offset,
        numcols*x_panel + numrows*y_panel + (numrows-1)*y_offset + (numcols-1)*x_offset,
    ]
    xs = [c[1] for c ∈ corners]
    ys = [c[2] for c ∈ corners]
    x1, x2 = minimum(xs), maximum(xs) 
    y1, y2 = minimum(ys), maximum(ys) 
    return abs((x2-x1)/(y2-y1))
end


function plot_chirality(Zs, sys;
    colorscheme=ColorSchemes.RdBu, clims = (-0.5, 0.5), offset_spacing = 1,
    numcols = nothing, texts = nothing, force_aspect = true, fig_kwargs...
)
    # Consolidate lattice info and panel layout
    numpanels = length(Zs)
    isnothing(numcols) && (numcols = numpanels)
    numrows = floor(Int, (numpanels - 1) / numcols ) + 1
    vecs = sys.lattice.lat_vecs
    v₁, v₂ = vecs[:,1], vecs[:,2]
    nx, ny = size(Zs[1])[1:2] 
    v₁, v₂ = Point3f(v₁), Point3f(v₂)
    x, y = [1., 0, 0], [0., 1, 0]
    x_offset = offset_spacing * (x ⋅ v₁)*x
    y_offset = -offset_spacing * (y ⋅ v₂)*y
    x_panel =  nx * v₁
    y_panel = -ny * v₂ 
    aspect = aspect_ratio(x_panel, y_panel, x_offset, y_offset, numrows, numcols)

    # Set up figure
    fig = Figure(; fig_kwargs...)
    if force_aspect 
        ax = Axis(fig[1,1:length(Zs)]; aspect)
    else
        ax = Axis(fig[1,1:length(Zs)])
    end
    hidespines!(ax); hidedecorations!(ax)

    # Plot panels
    plaq1(p) = GLMakie.Polygon(Point2f.([p, p+v₁, p+v₂]))
    plaq2(p) = GLMakie.Polygon(Point2f.([p+v₁, p+v₁+v₂, p+v₂]))
    for (i, Z) ∈ enumerate(Zs)
        r, c = fldmod1(i, numcols)
        v₀ = (c-1) * (x_panel + x_offset) + (r-1) * (y_panel + y_offset)
        Χ = plaquette_map(berry, Z)
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
        if !isnothing(texts)
            text!(ax, v₀[1], v₀[2]; text=texts[i])
        end
    end
    fig, ax
end

function plot_chirality(sys;
    colorscheme=ColorSchemes.RdBu, clims = (-0.5, 0.5), offset = 1,
    numcols = nothing, texts = nothing, fig_kwargs...
)
    plot_chirality([sys._coherents], sys; colorscheme, clims, offset, numcols, texts, fig_kwargs...)
end

function plot_spins_color!(ax, Zs, sys;
    arrowcolor = :firebrick2, arrowlength=1.5, arrowsize=0.3,
    linewidth=0.15, kwargs...
)
    points = GLMakie.Point3f0.(vec(sys.lattice))
    dipoles = Sunny.expected_spin.(Zs)
    vecs = GLMakie.Vec3f0.(vec(dipoles))
    points .-= points[1]  # Align with plaquettes, which start at (0,0,0)
    lengths = norm.(vecs)
    lengths ./= maximum(lengths)
    color = [(arrowcolor, length) for length ∈ lengths] # Set transparency according to dipole magnitude

    GLMakie.arrows!(ax, points, vecs;
        color, linewidth, arrowsize, lengthscale=arrowlength, kwargs...,
    )

    fig
end



################################################################################
# Animation functions 
################################################################################
function get_colors(colorscheme, x, clims=(-0.5, 0.5))
    nx, ny = size(x)[2:3] # First index is plaquette index
    colors = Array{ColorTypes.RGB{Float64}, 1}(undef, nx*ny*2)
    count = 1
    for c ∈ 1:ny, r ∈ 1:nx
        colors[count] = get(colorscheme, x[1,r,c,1,1], clims)
        count += 1
        colors[count] = get(colorscheme, x[2,r,c,1,1], clims)
        count += 1
    end
    colors
end

function animate_chirality(Zs, sys;
    filename = "chirality_anim.mp4", colorscheme=ColorSchemes.RdBu,
    clims=(-0.5, 0.5), framerate=30, skip_interval=1, text = nothing,
)
    # Consolidate lattice information
    lat_vecs = sys.lattice.lat_vecs
    v₁, v₂ = lat_vecs[:,1], lat_vecs[:,2]
    nx, ny = size(Zs)[1:2]
    v₁ = Point3f(v₁)
    v₂ = Point3f(v₂)
    corners = [(0, 0), nx*v₁, ny*v₂, nx*v₁ + ny*v₂]
    xs = [c[1] for c ∈ corners]
    ys = [c[2] for c ∈ corners]
    x1, x2 = minimum(xs), maximum(xs)
    y1, y2 = minimum(ys), maximum(ys)
    aspect = abs((x2-x1)/(y2-y1))


    # Set up plot elements using observables
    plaq1(p) = GLMakie.Polygon(Point2f.([p, p+v₁, p+v₂]))
    plaq2(p) = GLMakie.Polygon(Point2f.([p+v₁, p+v₁+v₂, p+v₂]))
    idx = Observable(1)
    Χ = @lift(plaquette_map(berry, Zs[:,:,:,:,$idx]))
    colors = @lift(get_colors(colorscheme, $Χ, clims))
    pgons = GLMakie.Polygon[]
    for c ∈ 1:ny, r ∈ 1:nx
        base = (r-1)*v₁ + (c-1)*v₂
        push!(pgons, plaq1(base))
        push!(pgons, plaq2(base))
    end

    # Set up figure
    fig = Figure(; resolution=(600,400))
    ax = Axis(fig[1,1]; aspect)
    hidespines!(ax); hidedecorations!(ax)
    poly!(ax, pgons; color=colors)
    if !isnothing(text)
        label = @lift(text[$idx])
        text!(ax, 0, 80.0, 0; text=label)
    end

    # Animate
    numsteps = size(Zs)[end]
    GLMakie.record(fig, filename, 1:skip_interval:numsteps; framerate) do i
        idx[] = i
    end
    nothing
end
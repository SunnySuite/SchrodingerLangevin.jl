function fig3()
    # Set up  model
    dims = (100, 100, 1)
    h = 15.35 # Field
    D = 19 # On-site anisotropy
    rng = MersenneTwister(111)
    sys = su3_skyrmion_model(dims; h, D, rng)

    # Set up integrator
    kT = 0.00 # Final temperature
    λ = 0.1 # Damping coefficient
    integrator = LangevinHeunP(sys, kT, λ)

    # Determine times to sample
    τs = [4., 16, 256]
    Δt = 0.01

    # Calculate dynamics
    println("Calculating dynamics...")
    rand!(sys)
    frames = []
    for i in 1:length(τs)
        dur = i == 1 ? τs[1] : τs[i] - τs[i-1] 
        numsteps = round(Int, dur/Δt)
        for _ ∈ 1:numsteps
            evolve!(integrator, Δt)
        end
        push!(frames, copy(sys._coherents))
    end

    # Create and save plot
    println("Generating figure...")
    fig, _ = plot_chirality(frames, sys;
        resolution=(1800,600), offset_spacing=10
    )
    save("fig3.png", fig)

    return fig
end
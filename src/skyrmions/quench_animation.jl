function quench_animation(;
    dur=256.0, sample_period=80, filename="chirality_animation.mp4"
)
    ## Set up model
    dims = (100, 100, 1)
    h = 15.35 # Field
    D = 19 # On-site anisotropy
    rng = MersenneTwister(111)
    sys = su3_skyrmion_model(dims; h, D, rng)
    rand!(sys)

    ## Set up integrator
    kT = 0.00 # Final temperature
    λ = 0.1 # Damping coefficient
    integrator = LangevinHeunP(sys, kT, λ)

    ## Run simulation
    println("Calculating dynamics...")

    Δt = 0.01
    numsteps = round(Int, dur/Δt)
    numframes = floor(Int, numsteps/sample_period)

    Zs = zeros(Sunny.CVec{3}, size(sys._coherents)..., numframes)
    c = 1
    for i ∈ 1:numsteps
        evolve!(integrator, Δt)
        if mod1(i, sample_period) == 1 && c <= numframes
            Zs[:,:,:,:,c] .= sys._coherents
            c += 1
        end
    end

    ## Generate animation
    println("Generating animation...")
    animate_chirality(Zs, sys; filename)
end
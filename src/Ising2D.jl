module Ising2D

using Statistics

# observables, including energy
include("Observables2D.jl")
export Energy, Magnetization

# utility functions
include("Utils2d.jl")
export RandomLattice

# spin flip sampler
include("samplers/SpinFlip2D.jl")
export SpinFlipMCMC, SpinFlipMCMC!

end # module end

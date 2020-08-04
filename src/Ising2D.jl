module Ising2D

using StaticArrays
using Statistics

# observables, including energy
include("Observables2D.jl")
export Energy, Magnetization

# utility functions
include("Utils2d.jl")
export RandomLattice

# spin flip sampler with Metropolis rates
include("samplers/Metropolis2D.jl")
export Metropolis, Metropolis!

end # module end

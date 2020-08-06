module Ising2D

using StaticArrays
using Statistics
using StatsBase

# observables, including energy
include("Observables2D.jl")
export Energy, Magnetization

# utility functions
include("Utils2d.jl")
export RandomLattice, RandomMLattice, GetNeighbors!

# spin flip sampler with Metropolis rates
include("samplers/Metropolis2D.jl")
export Metropolis, Metropolis!

# Wolff sampler
include("samplers/Wolff2D.jl")
export Wolff, Wolff!

end # module end

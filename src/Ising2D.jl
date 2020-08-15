module Ising2D

using StaticArrays
using Statistics
using StatsBase

# abstract types
include("types.jl")

# utility functions
include("utils2d.jl")
export RandomLattice, RandomMLattice, IsingOptions

# abstract sampler
include("sample.jl")
export sample_trajectory!, sample_trajectory

# observables, including energy
include("Observables2D.jl")
export Energy, Magnetization

# spin flip sampler with Metropolis rates
include("samplers/Metropolis2D.jl")
export Metropolis

# Wolff sampler
include("samplers/Wolff2D.jl")
export Wolff

end # module end

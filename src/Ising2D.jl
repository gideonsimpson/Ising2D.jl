module Ising2D

using Statistics
using StatsBase
using Printf

# abstract types
include("types.jl")

# utility functions
include("utils2D.jl")
export RandomLattice, RandomMLattice, IsingOptions

# abstract sampler
include("sample.jl")
export sample_trajectory!, sample_trajectory, sample_observables

# observables, including energy
include("Observables2D.jl")
export Energy, Magnetization, Clusters

# spin flip sampler with Metropolis rates
include("samplers/Metropolis2D.jl")
export Metropolis

# Wolff sampler
include("samplers/Wolff2D.jl")
export Wolff

end # module end

# Example use of the Metropolis sampler

using Plots
using Random
using Printf
using Statistics

push!(LOAD_PATH,"../src/")
using Ising2D

# Set lattice size, temperature, and number of iterations
N = 10;
β = 0.5;
niters = 10^6;

# Generate an initial random lattice
Random.seed!(100);
x0 = RandomLattice(N);

# Sample
x_trajectory = Metropolis(x0, β, niters);

# Compute statistics
e_trajectory = Energy.(x_trajectory)/N^2;
m_trajectory = Magnetization.(x_trajectory)/N^2;
@printf("β = %g, %d iterations\n\n", β, niters)
@printf("Mean Energy (per spin): %g, Std. Dev %g\n", mean(e_trajectory), std(e_trajectory))
@printf("Mean Magnetization (per spin): %g, Std. Dev %g\n", mean(m_trajectory), std(m_trajectory))

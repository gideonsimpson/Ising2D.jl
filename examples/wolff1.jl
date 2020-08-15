# Example use of the Wolff sampler

using Random
using Printf
using Statistics

push!(LOAD_PATH,"../src/")
using Ising2D

# Set lattice size, temperature, and number of iterations
N = 10;
β = 0.5;
J = 1.0; # lattice coupling
niters = 10^6;

# Generate an initial random lattice
Random.seed!(100);
x₀ = RandomLattice(L);
# Define the sampler
sampler = Wolff(β, J, L);
# Set sampler options
opts = IsingOptions(n_iters=niters,n_save_iters=1);

# Sample
x_trajectory = sample_trajectory(x₀, sampler, options=opts);

# Compute statistics
e_trajectory = Energy.(x_trajectory)/L^2;
m_trajectory = Magnetization.(x_trajectory)/L^2;
@printf("β = %g, %d iterations\n\n", β, niters)
@printf("Mean Energy (per spin): %g, Std. Dev %g\n", mean(e_trajectory), std(e_trajectory))
@printf("Mean Magnetization (per spin): %g, Std. Dev %g\n", mean(m_trajectory), std(m_trajectory))

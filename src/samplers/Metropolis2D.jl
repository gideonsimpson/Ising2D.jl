
struct Metropolis{TF<:AbstractFloat, TI<:Integer} <:AbstractSampler
    β::TF
    J::TF
    L::TI
end

"""
    `Metropolis(β, J, L)` - Set up the Metropolis sampler

### Fields
* `β` - Inverse temperature
* `J` - Coupling constant>0    
* `L` - Lattice size, L×L
"""
function Metropolis(β::TF, J::TF, L::TI) where{TF<:AbstractFloat, TI<:Integer}
    return Metropolis{TF, TI}(β, J, L)
end

mutable struct MetropolisState{TI<:Integer, Tx<:AbstractArray{TI}} <:AbstractSamplerState
    x::Tx
    x_proposal::Tx
end


function InitState!(x, sampler::Metropolis)

    return MetropolisState(x, deepcopy(x));
end

function InitState(x₀, sampler::Metropolis)

    return MetropolisState(deepcopy(x₀), deepcopy(x₀));
end


function UpdateState!(state::MetropolisState, sampler::Metropolis)

    i = rand(1:sampler.L);
    j = rand(1:sampler.L);
    state.x_proposal[i,j] *= -1;

    e = 0.0;

    k,l = mod1(i+1,sampler.L), j;
    e += - 0.5 * sampler.J * state.x_proposal[i,j] * state.x_proposal[k,l];
    k,l = mod1(i-1,sampler.L), j;
    e += - 0.5 * sampler.J * state.x_proposal[i,j] * state.x_proposal[k,l];
    k,l = i, mod1(j+1,sampler.L);
    e += - 0.5 * sampler.J * state.x_proposal[i,j] * state.x_proposal[k,l];
    k,l = i, mod1(j-1,sampler.L);
    e += - 0.5 * sampler.J * state.x_proposal[i,j] * state.x_proposal[k,l];

    # double once for double counting of bonds and again for the sign change
    ΔE = 2 * 2 * e;

    a = min(1, exp(-sampler.β * ΔE))

    if rand() <a
        state.x[i,j] *= -1;
    else
        state.x_proposal[i,j] *= -1;
    end
    state
end



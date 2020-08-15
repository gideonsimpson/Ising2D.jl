struct Wolff{TF<:AbstractFloat, TI<:Integer} <:AbstractSampler
    β::TF
    J::TF
    L::TI
    p::TF
end

"""
    `Wolff(β, J, L)` - Set up the Wolff sampler

### Fields
* `β` - Inverse temperature
* `J` - Coupling constant>0    
* `L` - Lattice size, L×L
"""
function Wolff(β::TF, J::TF, L::TI) where{TF<:AbstractFloat, TI<:Integer}
    return Wolff{TF, TI}(β, J, L, 1 - exp(-2 * β * J))
end


mutable struct WolffState{TI<:Integer, Tx<:AbstractArray{TI}} <:AbstractSamplerState
    x::Tx
    D::TI
end

function InitState!(x, sampler::Wolff)

    return WolffState(x, zero(eltype(x)));
end

function InitState(x₀, sampler::Wolff)

    return WolffState(deepcopy(x₀),zero(eltype(x₀)));
end

function UpdateState!(state::WolffState, sampler::Wolff)

    i = rand(1:sampler.L);
    j = rand(1:sampler.L);
    state.D = state.x[i,j];
    state.x[i,j] *=-1;

    k,l = mod1(i+1,sampler.L), j;
    if (state.x[k,l]==state.D)
        if (rand() < sampler.p)
            state.x[k,l] *=-1;
            FlipNeighbors!(state.x, k, l, sampler.p, state.D, sampler.L)
        end
    end

    k,l = mod1(i-1,sampler.L), j;
    if (state.x[k,l]==state.D)
        if (rand() < sampler.p)
            state.x[k,l] *=-1;
            FlipNeighbors!(state.x, k, l, sampler.p, state.D, sampler.L)
        end
    end

    k,l = i, mod1(j+1,sampler.L);
    if (state.x[k,l]==state.D)
        if (rand() < sampler.p)
            state.x[k,l] *=-1;
            FlipNeighbors!(state.x, k, l, sampler.p, state.D, sampler.L)
        end
    end

    k,l = i, mod1(j-1,sampler.L);
    if (state.x[k,l]==state.D)
        if (rand() < sampler.p)
            state.x[k,l] *=-1;
            FlipNeighbors!(state.x, k, l, sampler.p, state.D, sampler.L)
        end
    end
    state
end

"""
`FlipNeighbors!` - Recursively flip neighbors in the Wolff algorithm for 2D
square lattice Ising

### Fields
* `x` - state of the lattice
* `i` - i-coordinate of current spin
* `j` - j-coordinate of current spin
* `p` - Flipping probability
* `D` - Flipping direction
* `L` - Lattice size, N×N
"""
function FlipNeighbors!(x, i, j, p, D, L)

    k,l = mod1(i+1,L), j;
    if (x[k,l]==D)
        if (rand() < p)
            x[k,l] *=-1;
            FlipNeighbors!(x, k, l, p, D, L)
        end
    end

    k,l = mod1(i-1,L), j;
    if (x[k,l]==D)
        if (rand() < p)
            x[k,l] *=-1;
            FlipNeighbors!(x, k, l, p, D, L)
        end
    end

    k,l = i, mod1(j+1,L);
    if (x[k,l]==D)
        if (rand() < p)
            x[k,l] *=-1;
            FlipNeighbors!(x, k, l, p, D, L)
        end
    end

    k,l = i, mod1(j-1,L);
    if (x[k,l]==D)
        if (rand() < p)
            x[k,l] *=-1;
            FlipNeighbors!(x, k, l, p, D, L)
        end
    end
    x
end

"""
    `sample_trajectory!(x, sampler; options=IsingOptions())`

In place applciation of the `sampler` to `x`.  Number of iterations are
set using the `options` argument.

### Fields

* `x`         - Starting position for sampler, modified in place
* `sampler`   - Desired sampler
* `options`   - Sampling options, including number of iteration

"""
function sample_trajectory!(x::Tx, sampler::S; options=IsingOptions()) where {Tx, S<:AbstractSampler}

    state = InitState!(x, sampler);
    for _ in 1:options.n_iters
        UpdateState!(state, sampler);
    end
    x
end

"""
    `sample_trajectory(x₀, sampler; options=IsingOptions())`

Run the `sampler` starting at `x₀`.  Number of iterations and interval between
saves are set using the `options` argument.  For Metropolis samplers, the
running acceptance rates are also resturned.


### Fields

* `x`         - Starting position for sampler, modified in place
* `sampler`   - Desired sampler
* `options`   - Sampling options, including number of iteration

"""
function sample_trajectory(x₀::Tx, sampler::S; options=IsingOptions()) where {Tx,  S<:AbstractSampler}

    state = InitState(x₀, sampler);

    # allocate memory for samples
    samples = Tx[similar(x₀) for i = 1:options.n_save];
    save_index = 1;
    for n = 1:options.n_iters
        UpdateState!(state, sampler);
        if(mod(n,options.n_save_iters)==0)
            @. samples[save_index] = deepcopy(state.x);
            save_index+=1;
        end
    end
    return samples
end
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
saves are set using the `options` argument.

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

"""
    `sample_observables(x₀, sampler, observables; options=IsingOptions())`

Run the `sampler` starting at `x₀`.  Number of iterations and interval between
saves are set using the `options` argument. Returns the evaluation of hte
functions stored in observables at the desired save times.  This does not return
the full trajectory.


### Fields

* `x`         - Starting position for sampler, modified in place
* `sampler`   - Desired sampler
* `observables` - An array of scalar valued observable functions
* `options`   - Sampling options, including number of iteration

"""
function sample_observables(x₀::Tx, sampler::S, observables::TOV; options=IsingOptions()) where {Tx,  S<:AbstractSampler,TO<:Function, TOV<:AbstractVector{TO}}

    state = InitState(x₀, sampler);

    n_obs = length(observables);
    # allocate memory for samples
    observable_samples = zeros(n_obs, options.n_save);
    save_index = 1;
    for n = 1:options.n_iters
        UpdateState!(state, sampler);
        if(mod(n,options.n_save_iters)==0)
            for k in 1:n_obs
                observable_samples[k,save_index] = observables[k](state.x);
            end
            save_index+=1;
        end
    end
    return observable_samples
end
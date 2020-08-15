# old versions of the code

"""
`Metropolis` - Sample Ising on a N×N square lattice with spin flip proposals
and Metropolis rates

### Fields
* `x0` - Initial of the lattice
* `β` - Inverse temperature
* `niters` - Number of iterations to run
### Optional Fields
* `J = 1.0` - Coupling constant>0
* `nsave_iters = 1` - Number of iterations between saves
"""
function Metropolis(x0, β, niters; J=1.0, nsave_iters=1)

    N = size(x0)[1];
    x = copy(x0);
    xp = copy(x0);

    x_trajectory = typeof(x0)[];

    E = Energy(x, J = J);

    𝒩 = [zeros(MVector{2,Int8}) for _ in 1:4];

    for n in 1:niters
        i = rand(1:N);
        j = rand(1:N);
        xp[i,j] *= -1;

        # compute local energy contribution
        e = 0.0;
        GetNeighbors!(𝒩, i, j, N);
        for (k,l) in 𝒩
            e += - 0.5 * J * xp[i,j] * xp[k,l];
        end

        # double once for double counting of bonds and again for the sign change
        ΔE = 2 * 2 * e;

        a = min(1, exp(-β * ΔE))
        ζ = rand();
        if ζ < a
            x[i,j] *= -1;
            E += ΔE;
        else
            # restore state if rejected
            xp[i,j] *= -1;
        end
        if(mod(n,nsave_iters) == 0)
            push!(x_trajectory,copy(x));
        end
    end

    return x_trajectory
end

"""
`Metropolis!` - Sample Ising on a N×N square lattice with spin flip proposals
and Metropolis rates.  This is the in place version.

### Fields
* `x` - State of the lattice
* `β` - Inverse temperature
* `niters` - Number of iterations to run
### Optional Fields
* 'J=1.0' - Coupling constant>0
"""
function Metropolis!(x, β, niters; J=1.0)

    N = size(x)[1];
    xp = copy(x)

    E = Energy(x, J = J);

    𝒩 = [zeros(MVector{2,Int8}) for _ in 1:4];

    for n in 1:niters
        i = rand(1:N);
        j = rand(1:N);
        xp[i,j] *= -1;

        # compute local energy contribution
        e = 0.0;
        GetNeighbors!(𝒩, i, j, N);
        for (k,l) in 𝒩
            e += - 0.5 * J * xp[i,j] * xp[k,l];
        end

        # double once for double counting of bonds and again for the sign change
        ΔE = 2 * 2 * e;

        a = min(1, exp(-β * ΔE))
        ζ = rand();
        if ζ < a
            x[i,j] *= -1;
            E += ΔE;
        else
            # restore state if rejected
            xp[i,j] *= -1;
        end
    end
    x
end


"""
`Wolff` - Sample Ising on a N×N square lattice with the Wolff cluster algorithm

### Fields
* `x0` - Initial state of the lattice
* `β` - Inverse temperature
* `niters` - Number of iterations to run
### Optional Fields
* 'J = 1.0' - Coupling constant>0
* `nsave_iters = 1` - Number of iterations between saves
"""
function Wolff(x0, β, niters; J = 1.0, nsave_iters=1)

    p = 1 - exp(-2 * β * J);

    N = size(x0)[1];
    x = copy(x0);

    x_trajectory = typeof(x0)[];

    for n in 1:niters
        i = rand(1:N);
        j = rand(1:N);

        D = x[i,j];
        x[i,j] *=-1;

        k,l = mod1(i+1,N), j;
        if (x[k,l]==D)
            if (rand() < p)
                x[k,l] *=-1;
                FlipNeighbors!(x, k, l, p, D, N)
            end
        end

        k,l = mod1(i-1,N), j;
        if (x[k,l]==D)
            if (rand() < p)
                x[k,l] *=-1;
                FlipNeighbors!(x, k, l, p, D, N)
            end
        end

        k,l = i, mod1(j+1,N);
        if (x[k,l]==D)
            if (rand() < p)
                x[k,l] *=-1;
                FlipNeighbors!(x, k, l, p, D, N)
            end
        end

        k,l = i, mod1(j-1,N);
        if (x[k,l]==D)
            if (rand() < p)
                x[k,l] *=-1;
                FlipNeighbors!(x, k, l, p, D, N)
            end
        end
        if(mod(n,nsave_iters) == 0)
            push!(x_trajectory,copy(x));
        end
    end

    return x_trajectory
end

"""
`Wolff!` - Sample Ising on a N×N square lattice with the Wolff cluster
algorithm.  In place version.

### Fields
* `x` - State of the lattice
* `β` - Inverse temperature
* `niters` - Number of iterations to run
### Optional Fields
* 'J=1.0' - Coupling constant>0
"""
function Wolff!(x, β, niters; J = 1.0)

    p = 1 - exp(-2 * β * J);

    N = size(x)[1];

    for n in 1:niters
        i = rand(1:N);
        j = rand(1:N);

        D = x[i,j];
        x[i,j] *=-1;

        k,l = mod1(i+1,N), j;
        if (x[k,l]==D)
            if (rand() < p)
                x[k,l] *=-1;
                FlipNeighbors!(x, k, l, p, D, N)
            end
        end

        k,l = mod1(i-1,N), j;
        if (x[k,l]==D)
            if (rand() < p)
                x[k,l] *=-1;
                FlipNeighbors!(x, k, l, p, D, N)
            end
        end

        k,l = i, mod1(j+1,N);
        if (x[k,l]==D)
            if (rand() < p)
                x[k,l] *=-1;
                FlipNeighbors!(x, k, l, p, D, N)
            end
        end

        k,l = i, mod1(j-1,N);
        if (x[k,l]==D)
            if (rand() < p)
                x[k,l] *=-1;
                FlipNeighbors!(x, k, l, p, D, N)
            end
        end

    end

    x
end

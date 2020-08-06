"""
`Wolff` - Sample Ising on a N×N square lattice with the Wolff cluster algorithm

### Fields
* `x0` - Initial state of the lattice
* `β` - Inverse temperature
* `niters` - Number of iterations to run
### Optional Fields
* 'J=1.0' - Coupling constant>0
"""
function Wolff(x0, β, niters; J = 1.0)

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

        push!(x_trajectory,copy(x));
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

"""
`FlipNeighbors!` - Recursively flip neighbors in the Wolff algorithm for 2D
square lattice Ising

### Fields
* `x` - state of the lattice
* `i` - i-coordinate of current spin
* `j` - j-coordinate of current spin
* `p` - Flipping probability
* `D` - Flipping direction
* `N` - Lattice size, N×N
"""
function FlipNeighbors!(x, i, j, p, D, N)

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

    x
end

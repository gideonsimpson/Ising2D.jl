
"""
`Metropolis` - Sample Ising on a N×N square lattice with spin flip proposals
and Metropolis rates

### Fields
* `x0` - Initial of the lattice
* `β` - Inverse temperature
* `niters` - Number of iterations to run
### Optional Fields
* 'J=1.0' - Coupling constant>0
"""
function Metropolis(x0, β, niters; J=1.0)

    N = size(x0)[1];
    x = copy(x0);
    xp = copy(x0);

    x_trajectory = typeof(x0)[];

    E = Energy(x, J = J);

    for n in 1:niters
        i = rand(1:N);
        j = rand(1:N);
        xp[i,j] *= -1;

        # compute local energy contribution
        e = 0.0;

        k = mod1(i+1,N);
        l = j;
        e += - 0.5 * J * xp[i,j] * xp[k,l];

        k = mod1(i-1,N);
        l = j;
        e += - 0.5 * J * xp[i,j] * xp[k,l];

        k = i
        l = mod1(j+1,N);;
        e += - 0.5 * J * xp[i,j] * xp[k,l];

        k = i
        l = mod1(j-1,N);
        e += - 0.5 * J * xp[i,j] * xp[k,l];

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
        push!(x_trajectory,copy(x));
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

    for n in 1:niters
        i = rand(1:N);
        j = rand(1:N);
        xp[i,j] *= -1;

        # compute local energy contribution
        e = 0.0;

        k = mod1(i+1,N);
        l = j;
        e += - 0.5 * J * xp[i,j] * xp[k,l];

        k = mod1(i-1,N);
        l = j;
        e += - 0.5 * J * xp[i,j] * xp[k,l];

        k = i
        l = mod1(j+1,N);;
        e += - 0.5 * J * xp[i,j] * xp[k,l];

        k = i
        l = mod1(j-1,N);
        e += - 0.5 * J * xp[i,j] * xp[k,l];

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

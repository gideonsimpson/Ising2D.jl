
"""
`Energy` - Compute the energy

E = -0.5 * J ∑ xᵢxⱼ

over a periodic N×N square lattice.  Note that this counts each bond twice,
hence the factor of 0.5

### Fields
* `x` - State of the lattice
### Optional Fields
* 'J=1' - Coupling constant>0
"""
function Energy(x; J=1.0)
    E = 0.0;

    N= size(x)[1];
    for i in 1:N, j in 1:N
        # sum of energy contributions

        k = mod1(i+1,N);
        l = j;
        E += - 0.5 * J * x[i,j] * x[k,l];

        k = mod1(i-1,N);
        l = j;
        E += - 0.5 * J * x[i,j] * x[k,l];

        k = i
        l = mod1(j+1,N);;
        E += - 0.5 * J * x[i,j] * x[k,l];

        k = i
        l = mod1(j-1,N);
        E += - 0.5 * J * x[i,j] * x[k,l];
    end
    return E
end


"""
`Magnetization` - Compute the magnetization

M = ∑ xᵢ

over a periodic N×N square lattice.

### Fields
* `x` - State of the lattice
"""
function Magnetization(x)
    M = sum(x);
    return M
end

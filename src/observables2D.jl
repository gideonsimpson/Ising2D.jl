
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
function Energy(x::Tx; J::TF=1.0) where {TI<:Integer, Tx<:AbstractArray{TI}, TF<:AbstractFloat}
    E = 0.0;

    L= size(x)[1];
    for i in 1:L, j in 1:L
        # sum of energy contributions
        k, l = mod1(i+1,L), j;
        E += - 0.5 * J * x[i,j] * x[k,l];

        k,l = mod1(i-1,L), j;
        E += - 0.5 * J * x[i,j] * x[k,l];

        k, l = i, mod1(j+1,L);;
        E += - 0.5 * J * x[i,j] * x[k,l];

        k, l = i, mod1(j-1,L);
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
function Magnetization(x::Tx) where {TI<:Integer, Tx<:AbstractArray{TI}}
    return sum(x)
end


"""
`ClusterNeighbors!` - Mark on the cluster map, `c`, where the current cluster
continues.


### Fields
* `c` - Cluster map
* `x` - State of the lattice
* `i,j` - Indices of the current site
* `D` - Direction of the current cluster (+/-1)
* `id` - Cluster id
* `L` - Lattice size, L × L
"""
function ClusterNeighbors!(c, x, i, j, D, id, L)

    k,l = mod1(i+1,L), j;
    if (x[k,l] == D && c[k,l]==0)
        c[k,l] = id;
        ClusterNeighbors!(c, x, k, l, D, id, L)
    end

    k,l = mod1(i-1,L), j;
    if (x[k,l]==D && c[k,l]==0)
        c[k,l] = id;
        ClusterNeighbors!(c, x, k, l, D, id, L)
    end

    k,l = i, mod1(j+1,L);
    if (x[k,l]==D && c[k,l]==0)
        c[k,l] = id;
        ClusterNeighbors!(c, x, k, l, D, id, L)
    end

    k,l = i, mod1(j-1,L);
    if (x[k,l]==D && c[k,l]==0)
        c[k,l] = id;
        ClusterNeighbors!(c, x, k, l, D, id, L)
    end
    c
end


"""
`Cluster` - Find all the clusters in the lattice.  Note that the numbering is
arbitrary.


### Fields
* `x` - State of the lattice
"""
function Clusters(x::Tx) where {TI<:Integer, Tx<:AbstractArray{TI}}
    L= size(x)[1];
    id = 0;
    c = zeros(Int, L, L);

    for i in 1:L, j in 1:L
        if(c[i,j]==0)
            id+=1;
            c[i,j] = id;
            D = x[i,j];
            ClusterNeighbors!(c, x, i, j, D, id, L)
        end
    end

    return c
end


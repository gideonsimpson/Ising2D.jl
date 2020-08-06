
"""
`RandomLattice` - Generate a random square lattice

### Fields
* `N` - Lattice size, N×N
"""
function RandomLattice(N)
    return Int8.(2*rand(Bool, N, N).-1)
end

"""
`RandomMLattice` - Generate a random square lattice with magnetization M.  This
assumes that M is compatible with the allowable magenetizations of the lattice

### Fields
* `N` - Lattice size, N×N
* `M` - Target magnetization, M
"""
function RandomMLattice(N,M)
    
    if M ≥ 0
        σ = Int8(1);
    else
        σ = Int8(-1);
    end
    x = σ  * ones(Int8, N, N);

    # compute number of sites to flip
    n_flip = (N^2 -abs(M)) ÷ 2;
    # sample the sites
    flip_sites = sample(collect(Iterators.product(1:N,1:N))[:],n_flip,replace=false)
    # flip
    for site in flip_sites
        x[site[1],site[2]] *=Int8(-1);
    end

    return x

end

"""
`GetNeighbors!` - Generate a list of neighbors

### Fields
* `𝒩` - Neighbor set to be populated.  Must be a mutable type
* `i` - i-index
* `j` - j-index
* `N` - Lattice size
"""
function GetNeighbors!(𝒩, i, j, N)
    𝒩[1][1] = mod1(i+1,N);
    𝒩[1][2] = j;

    𝒩[2][1] = mod1(i-1,N);
    𝒩[2][2] = j;

    𝒩[3][1] = i;
    𝒩[3][2] = mod1(j+1,N);

    𝒩[4][1] = i;
    𝒩[4][2] = mod1(j-1,N);

    𝒩
end

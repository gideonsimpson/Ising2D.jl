
"""
`RandomLattice` - Generate a random square lattice

### Fields
* `N` - Lattice size, N×N
"""
function RandomLattice(N)
    return Int8.(2*rand(Bool, N, N).-1)
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

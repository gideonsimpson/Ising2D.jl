
"""
`RandomLattice` - Generate a random square lattice

### Fields
* `N` - Lattice size, N×N
"""
function RandomLattice(N)
    return Int8.(2*rand(Bool, N, N).-1)
end

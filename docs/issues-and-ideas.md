- How do we connect the types from MeshCore to the Elfel element types?
- TODO: Do I need an abstract type for the finite element?
- Why does the FE iterator allocate? 
p = _storedofs!(it._dofs, p, state, it._irs[i], it._flds[i])
Is it the abstract types of the IRs and fields?
- For multiple FE spaces, one could have an iterator which simultaneously iterates multiple spaces (invoking their own iterators).
```
struct It1
    a
end
_update!(it::It1, state) = (it.a[state] = state)
function Base.iterate(it::It1, state = 1)
    if state > length(it.a)
        return nothing
    else
        return (_update!(it, state), state+1)
    end
end
Base.length(it::It1)  = length(it.a)
struct It2
    b
end
_update!(it::It2, state) = (it.b[state] = -state)
function Base.iterate(it::It2, state = 1)
    if state > length(it.b)
        return nothing
    else
        return (_update!(it, state), state+1)
    end
end
Base.length(it::It2)  = length(it.b)
struct It
    it1::It1
    it2::It2
end
function _update!(it::It, state) 
    iterate(it.it1, state)
    iterate(it.it2, state)
    it
end
function Base.iterate(it::It, state = 1)
    if state > length(it.it1)
        return nothing
    else
        return (_update!(it, state), state+1)
    end
end
Base.length(it::It2)  = length(it.b)

N = 3
it = It(It1(fill(0, N)), It2(fill(0, N)))
for i in it
    @show it.it1
    @show it.it2
end
```

- FE space could figure out basis function classification: which entity is it attached to, which one is it in serial order, and anything else, and share it with iterators (especially the quadrature-point iterator). Done.

- Some tests with strain-displacement matrices.
```
module m
using BenchmarkTools
using StaticArrays
using LinearAlgebra

B = (g -> SVector{3}((g[1], 0, g[2])),
   g -> SVector{3}((0, g[2], g[1])))
function Bf(g, k)
    if k == 1
        SVector{3}((g[1], 0, g[2]))
    else
        SVector{3}((0, g[2], g[1]))
    end
end
c = [1, 2, 1, 2, 1, 2]
g = SVector{2}(rand(2))

# With allocations
@btime begin 
    v = 0.0
    for i in 1:6
        Bi = $B[$c[i]]($g)
        for j in 1:6
            Bj = $B[$c[j]]($g)
            v += dot(Bi, Bj)
        end
    end
    v
end

# No allocations
@btime begin 
    v = 0.0
    ci = 1
    for i in 1:6
        Bi = $B[ci]($g)
        ci = ci == 1 ? 2 : 1 
        cj = 1
        for j in 1:6
            Bj = $B[cj]($g)
            v += dot(Bi, Bj)
            cj = cj == 1 ? 2 : 1 
        end
    end
    v
end

# No allocations
@btime begin 
    v = 0.0
    for i in 1:6
        Bi = $Bf($g, $c[i])
        for j in 1:6
            Bj = $Bf($g, $c[j])
            v += dot(Bi, Bj)
        end
    end
    v
end

function test()
# With allocations
begin 
    v = 0.0
    for i in 1:6
        Bi = B[c[i]](g)
        for j in 1:6
            Bj = B[c[j]](g)
            v += dot(Bi, Bj)
        end
    end
    @show v
end

# No allocations
begin 
    v = 0.0
    ci = 1
    for i in 1:6
        Bi = B[ci](g)
        ci = ci == 1 ? 2 : 1 
        cj = 1
        for j in 1:6
            Bj = B[cj](g)
            v += dot(Bi, Bj)
            cj = cj == 1 ? 2 : 1 
        end
    end
    @show v
end
end
test()

end
```

- L2 element may be incapable of representing its geometry because it does not necessarily have nodal basis functions. What do we do about that?



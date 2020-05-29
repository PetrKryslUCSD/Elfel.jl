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
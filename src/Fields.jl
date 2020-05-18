module Fields

import Base: setindex!
using StaticArrays

struct FEFieldVector{VT}
    _v::Vector{VT}
end

function (o::FEFieldVector{VT})(j::Int64) where {VT} 
    o._v[j]
end

function setindex!(o::FEFieldVector{VT}, X, inds...) where {VT}
    o._v[inds...] = X
    return o
end

"""
    FEField{N, T}

Structure to store degree of freedom numbers and values.
The vectors are addressed by the entity number.
"""
mutable struct FEField{N, T}
    dofvecs::FEFieldVector{SVector{N, T}}
    isknown::FEFieldVector{SVector{N, Bool}}
    dofnums::FEFieldVector{SVector{N, Int64}}
    nunknowns::Int64
end

function FEField{N, T}(nent::Int64) where {N, T}
    z = SVector{N, T}([zero(T) for i in 1:N])
    dofvecs = FEFieldVector([z for idx in 1:nent]) 
    z = SVector{N, Bool}([false for i in 1:N])
    isknown = FEFieldVector([z for idx in 1:nent]) 
    z = SVector{N, Int64}([0 for i in 1:N])
    dofnums = FEFieldVector([z for idx in 1:nent]) 
    return FEField{N, T}(dofvecs, isknown, nums, 0)
end

"""
    ndofsperentity(self::FEField{N, T}) where {N, T}

How many degrees of freedom per entity?
"""
ndofsperentity(self::FEField{N, T}) where {N, T} = N

"""
    nentities(self::FEField{N, T}) where {N, T} 

Number of entities in the field.
"""
nentities(self::FEField{N, T}) where {N, T} = length(self.dofvecs._v)

ndofs(self::FEField{N, T}) where {N, T} = nentities(self) * ndofsperentity(self)

"""
    numberdofs!(self::AbstractField)

Number the degrees of freedom.

The free components in the field are numbered consecutively. No effort is
made to optimize the numbering in any way. If you'd like to optimize the
numbering of the degrees of freedom, use the above form that sets the
permutation of the degrees of freedom, or the permutation of the nodes.
"""
function numberdofs!(self::FEField{N, T}) where {N, T}
    num = zero(typeof(self.nunknowns))
    nents = nentities(self)
    ndpe = ndofsperentity(self)
    for i in 1:nents
        n = MVector(self.nums(i))
        for j in 1:ndpe
            if !self.isknown._v[i][j] # free degree of freedom
                num = num + 1
                n[j] = num
            end
        end
        self.nums[i] = n
    end
    self.nunknowns = num
    for i in 1:nents
        n = MVector(self.nums(i))
        for j in 1:ndpe
            if self.isknown._v[i][j] # free degree of freedom
                num = num + 1
                n[j] = num
            end
        end
        self.nums[i] = n
    end
    return  self
end

"""
    setebc!(self::FEField{N, T}, eid, comp, val::T) where {N, T}

Set the EBCs (essential boundary conditions).

- `eid`  = entity identifier,
- `comp` = integer, which  degree of freedom (component),
- `val`  = value of type T
"""
function setebc!(self::FEField{N, T}, eid, comp, val::T) where {N, T}
    ik = MVector(self.isknown(eid))
    ik[comp] = true
    self.isknown[eid] = ik
    d = MVector(self.dofvecs(eid))
    d[comp] = val
    self.dofvecs[eid] = d
    return  self
end

"""
    gathersysvec(self::FEField{N, T}) where {N, T}

Gather values from the field for the whole system vector.
"""
function gathersysvec(self::FEField{N, T}) where {N, T}
    nents = nentities(self)
    ndpe = ndofsperentity(self)
    vec = fill(zero(T), ndofs(self))
    for i in 1:nents
        en = self.dofnums(i)
        for j in 1:ndpe
            vec[en[j]] = self.dofvecs(i)[j]
        end
    end
    return vec
end

end

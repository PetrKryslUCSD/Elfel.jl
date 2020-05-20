module FEFields

using StaticArrays
using MeshCore: nshapes

"""
    FEField{N, T, S, IT}

Type of FE field.
"""
mutable struct FEField{N, T, S, IT}
    shapecollection::S
    dofnums::Vector{SVector{N, IT}}
    isdatum::Vector{SVector{N, Bool}}
    dofvals::Vector{SVector{N, T}}
    nunknowns::Int64
    function FEField(::Val{N}, ::Type{T}, ::Type{IT}, sc::S) where {N, T, S, IT}
        nv = nshapes(sc)
        z = fill(zero(IT), N)
        dofnums = [SVector{N}(z) for i in 1:nv]
        z = fill(false, N)
        isdatum = [z for i in 1:nv]
        z = fill(zero(T), N)
        dofvals = [SVector{N}(z) for i in 1:nv]
        return new{N, T, S, IT}(sc, dofnums, isdatum, dofvals)
    end
end

nterms(fef::FEField) = nshapes(fef.shapecollection)
ndofsperterm(fef::FEField{N}) where {N} = N
ndofs(fef::FEField) = nterms(fef) * ndofsperterm(fef)

"""
    numberdofs!(self::FEField)

Number the degrees of freedom.

The unknown degrees of freedom in the field are numbered consecutively. Then the
datum degrees of freedom are numbered, again consecutively. 

No effort is made to optimize the numbering in any way. 
"""
function numberdofs!(self::FEField) 
    num = zero(typeof(self.nunknowns))
    nt = nterms(self)
    ndpt = ndofsperterm(self)
    for i in 1:nt
        n = MVector(self.dofnums[i])
        for j in 1:ndpt
            if !self.isdatum[i][j] # unknown degree of freedom
                num = num + 1
                n[j] = num
            end
        end
        self.dofnums[i] = n
    end
    self.nunknowns = num
    for i in 1:nt
        n = MVector(self.dofnums[i])
        for j in 1:ndpt
            if self.isdatum[i][j] # datum degree of freedom
                num = num + 1
                n[j] = num
            end
        end
        self.dofnums[i] = n
    end
    return  self
end

"""
    setebc!(self::FEField{N, T}, eid, comp, val::T) where {N, T}

Set the EBCs (essential boundary conditions).

- `tid`  = term identifier (serial number),
- `comp` = which  degree of freedom in the term,
- `val`  = value of type T
"""
function setebc!(self::FEField, tid, comp, val::T) where {T}
    ik = MVector(self.isdatum[tid])
    ik[comp] = true
    self.isdatum[tid] = ik
    d = MVector(self.dofvals[tid])
    d[comp] = val
    self.dofvals[tid] = d
    return  self
end

"""
    gathersysvec(self::FEField{N, T}) where {N, T}

Gather values from the field for the whole system vector.
"""
function gathersysvec(self::FEField{N, T}) where {N, T}
    nt = nterms(self)
    ndpt = ndofsperterm(self)
    vec = fill(zero(T), nt * ndpt)
    for i in 1:nt
        en = self.dofnums[i]
        for j in 1:ndpt
            vec[en[j]] = self.dofvals[i][j]
        end
    end
    return vec
end

end

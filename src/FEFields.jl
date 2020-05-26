module FEFields

using StaticArrays
using ..FElements: nfeatofdim, ndofsperfeat
import ..FElements: ndofsperelem

mutable struct FEField{N, T, IT}
    dofnums::Vector{SVector{N, IT}}
    isdatum::Vector{SVector{N, Bool}}
    dofvals::Vector{SVector{N, T}}
    nunknowns::Int64

    function FEField(::Val{N}, ::Type{T}, ::Type{IT}, ns) where {N, T, IT}
        z = fill(zero(IT), N)
        dofnums = [SVector{N}(z) for i in 1:ns]
        z = fill(false, N)
        isdatum = [z for i in 1:ns]
        z = fill(zero(T), N)
        dofvals = [SVector{N}(z) for i in 1:ns]
        return new{N, T, IT}(dofnums, isdatum, dofvals, 0)
    end
end

doftype(fef::FEField{N, T, IT}) where {N, T, IT} = T
nterms(fef::FEField) = length(fef.dofnums)
ndofsperterm(fef::FEField{N}) where {N} = N
ndofs(fef::FEField) = nterms(fef) * ndofsperterm(fef)
nunknowns(fef::FEField) = fef.nunknowns

function setebc!(self::FEField, tid, comp, val::T) where {T}
    ik = MVector(self.isdatum[tid])
    ik[comp] = true
    self.isdatum[tid] = ik
    d = MVector(self.dofvals[tid])
    d[comp] = val
    self.dofvals[tid] = d
    return  self
end

function numberdofs!(f::FEField) 
    num = zero(typeof(f.nunknowns))
    nt = nterms(f)
    ndpt = ndofsperterm(f)
    for i in 1:nt
        n = MVector(f.dofnums[i])
        for j in 1:ndpt
            if (!f.isdatum[i][j]) # unknown degree of freedom
                num = num + 1
                n[j] = num
            end
        end
        f.dofnums[i] = n
    end
    f.nunknowns = num
    for i in 1:nt
        n = MVector(f.dofnums[i])
        for j in 1:ndpt
            if f.isdatum[i][j] # datum degree of freedom
                num = num + 1
                n[j] = num
            end
        end
        f.dofnums[i] = n
    end
    return  f
end

function gathersysvec!(vec, self::FEField)
    nt = nterms(self)
    ndpt = ndofsperterm(self)
    # vec = fill(zero(eltype(self.dofvals[1])), nt * ndpt)
    for i in 1:nt
        en = self.dofnums[i]
        for j in 1:ndpt
            vec[en[j]] = self.dofvals[i][j]
        end
    end
    return vec
end

function scattersysvec!(self::FEField, v)
    nt = nterms(self)
    ndpt = ndofsperterm(self)
    for i in 1:nt
        en = self.dofnums[i]
        for j in 1:ndpt
            self.dofvals[i][j] = v[en[j]]
        end
    end
    return self
end

end

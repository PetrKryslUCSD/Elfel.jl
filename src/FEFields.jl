module FEFields

using StaticArrays

mutable struct FEField{N, T, IT}
    dofnums::Vector{SVector{N, IT}}
    isdatum::Vector{SVector{N, Bool}}
    dofvals::Vector{SVector{N, T}}

    function FEField(::Val{N}, ::Type{T}, ::Type{IT}, ns) where {N, T, IT}
        z = fill(zero(IT), N)
        dofnums = [SVector{N}(z) for i in 1:ns]
        z = fill(false, N)
        isdatum = [z for i in 1:ns]
        z = fill(zero(T), N)
        dofvals = [SVector{N}(z) for i in 1:ns]
        return new{N, T, IT}(dofnums, isdatum, dofvals)
    end
end

doftype(fef::FEField{N, T, IT}) where {N, T, IT} = T
dofnumtype(fef::FEField{N, T, IT}) where {N, T, IT} = IT
nterms(fef::FEField) = length(fef.dofnums)
ndofsperterm(fef::FEField{N}) where {N} = N
ndofs(fef::FEField) = nterms(fef) * ndofsperterm(fef)

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
    numberfreedofs!(f::FEField, firstnum = 1) 

Number the unknowns in the field, starting from the one supplied on input.

Note: The data degrees of freedom have their numbers zeroed out.
"""
function numberfreedofs!(f::FEField, firstnum = 1) 
    lnum = fnum = zero(dofnumtype(f)) + firstnum
    tnum = 0
    for i in 1:nterms(f)
        n = MVector(f.dofnums[i])
        for j in 1:ndofsperterm(f)
            if (!f.isdatum[i][j]) # unknown degree of freedom
                n[j] = lnum
                tnum = tnum + 1
                lnum = lnum + 1
            else
                n[j] = zero(dofnumtype(f))
            end
        end
        f.dofnums[i] = n
    end
    return  f
end

"""
    numberdatadofs!(f::FEField, firstnum = 1)

Number the data degrees of freedom in the field. Start from the
number supplied on input.

Note: The free degrees of freedom are numbered first.
"""
function numberdatadofs!(f::FEField, firstnum = 1) 
    num = zero(dofnumtype(f)) + firstnum
    for i in 1:nterms(f)
        n = MVector(f.dofnums[i])
        for j in 1:ndofsperterm(f)
            if (f.isdatum[i][j]) # known (data) degree of freedom
                n[j] = num
                num = num + 1
            end
        end
        f.dofnums[i] = n
    end
    return  f
end

function freedofnums(f::FEField) 
    tnum = lnum = zero(dofnumtype(f))
    fnum = typemax(dofnumtype(f))
    for i in 1:nterms(f)
        for j in 1:ndofsperterm(f)
            # unknown degree of freedom
            if (!f.isdatum[i][j]) 
                tnum = tnum + 1
                lnum = max(lnum, f.dofnums[i][j])
                fnum = min(fnum, f.dofnums[i][j])
            end
        end
    end
    return  (fnum, lnum, tnum)
end

highestfreedofnum(f::FEField) = freedofnums(f)[2]

function datadofnums(f::FEField) 
    tnum = lnum = zero(dofnumtype(f))
    fnum = typemax(dofnumtype(f))
    for i in 1:nterms(f)
        for j in 1:ndofsperterm(f)
            # known degree of freedom
            if (f.isdatum[i][j]) 
                tnum = tnum + 1
                lnum = max(lnum, f.dofnums[i][j])
                fnum = min(fnum, f.dofnums[i][j])
            end
        end
    end
    return  (fnum, lnum, tnum)
end

highestdatadofnum(f::FEField) = datadofnums(f)[2]

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
    vl = length(v)
    nt = nterms(self)
    ndpt = ndofsperterm(self)
    for i in 1:nt
        en = self.dofnums[i]
        for j in 1:ndpt
            if en[j] <= vl
                d = MVector(self.dofvals[i])
                d[j] = v[en[j]]
                self.dofvals[i] = d
            end
        end
    end
    return self
end

end

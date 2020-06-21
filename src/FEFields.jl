module FEFields

using StaticArrays

"""
    FEField{N, T, IT}

Type of a finite element field. Parameterized with
- `N`: number of degrees of freedom per entity, 
- `T`: type of the degree of freedom value, 
- `IT`: type of the index (integer value). This describes the serial numbers
  of the degrees of freedom.
"""
mutable struct FEField{N, T, IT}
    dofnums::Vector{SVector{N, IT}}
    isdatum::Vector{SVector{N, Bool}}
    dofvals::Vector{SVector{N, T}}

    """
        FEField(::Val{N}, ::Type{T}, ::Type{IT}, ntrms) where {N, T, IT}
    
    Construct a field with the given number of terms `ntrms`.
    """
    function FEField(::Val{N}, ::Type{T}, ::Type{IT}, ntrms) where {N, T, IT}
        z = fill(zero(IT), N)
        dofnums = [SVector{N}(z) for i in 1:ntrms]
        z = fill(false, N)
        isdatum = [z for i in 1:ntrms]
        z = fill(zero(T), N)
        dofvals = [SVector{N}(z) for i in 1:ntrms]
        return new{N, T, IT}(dofnums, isdatum, dofvals)
    end
end

"""
    doftype(fef::FEField{N, T, IT}) where {N, T, IT}

Type of a degree of freedom value.
"""
doftype(fef::FEField{N, T, IT}) where {N, T, IT} = T

"""
    dofnumtype(fef::FEField{N, T, IT}) where {N, T, IT}

Type of the index (serial number) of the degree of freedom. Integer.
"""
dofnumtype(fef::FEField{N, T, IT}) where {N, T, IT} = IT
"""
    nterms(fef::FEField)

Number of terms in the field.
"""
nterms(fef::FEField) = length(fef.dofnums)

"""
    ndofsperterm(fef::FEField{N}) where {N}

Number of degrees of freedom per term of the field.
"""
ndofsperterm(fef::FEField{N}) where {N} = N

"""
    ndofs(fef::FEField)

Total number of degrees of freedom in the field.
"""
ndofs(fef::FEField) = nterms(fef) * ndofsperterm(fef)

"""
    setebc!(self::FEField, tid, comp, val::T) where {T}

Set the value of one particular degree of freedom to a given number.

- `tid`: which term, 
- `comp`: which component of the term, 
- `val`: value to which the degree of freedom should be set.
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

Note: The free degrees of freedom must be numbered first.
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

"""
    freedofnums(f::FEField) 

Collect information about unknown (free) degree of freedom numbers.

First number, last number, and the total number of degrees of freedom are
returned as a tuple.
"""
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


"""
    highestfreedofnum(f::FEField)

Compute the highest serial number of a free degree of freedom in the field.
"""
highestfreedofnum(f::FEField) = freedofnums(f)[2]

"""
    datadofnums(f::FEField) 

Collect information about known (data) degree of freedom numbers.

First number, last number, and the total number of degrees of freedom are
returned as a tuple.
"""
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

"""
    highestdatadofnum(f::FEField)

Compute the highest serial number of a datum degree of freedom in the field.
"""
highestdatadofnum(f::FEField) = datadofnums(f)[2]

"""
    gathersysvec!(vec, self::FEField)

Gather system vector contributions from the field.
"""
function gathersysvec!(vec, self::FEField)
    nt = nterms(self)
    ndpt = ndofsperterm(self)
    for i in 1:nt
        en = self.dofnums[i]
        for j in 1:ndpt
            vec[en[j]] = self.dofvals[i][j]
        end
    end
    return vec
end

"""
    scattersysvec!(self::FEField, v)

Scatter a system vector into the field.
"""
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

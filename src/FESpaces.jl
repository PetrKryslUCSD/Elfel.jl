module FESpaces

using StaticArrays
using MeshCore
using MeshCore: nshapes, indextype, nrelations, retrieve, manifdim, IncRel
using MeshKeeper: Mesh, baseincrel, increl
using ..FElements: nfeatofdim, ndofsperfeat, ndofsperelem

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

struct FESpace{FET}
    fe::FET
    mesh::Mesh
    _irsfields::Dict

    function FESpace(::Type{T}, fe::FET, mesh) where {T, FET}
        baseir = baseincrel(mesh)
        _irsfields = _makefields(T, indextype(baseir), fe, mesh)
        return new{FET}(fe, mesh, _irsfields)
    end
end

function _makefields(::Type{T}, ::Type{IT}, fe, mesh) where {T, IT} 
    _irsfields= Dict()
    for m in 0:1:manifdim(fe.sd)
        if ndofsperfeat(fe, m) > 0
            fir = increl(mesh, (manifdim(fe.sd), m))
            fld = FEField(Val(ndofsperfeat(fe, m)), T, IT, nshapes(fir.right))
            _irsfields[m] = (fir, fld)
        end 
    end
    return _irsfields
end

#= TODO is it more natural to have access to the geometry from the font element space or from the iterator? =#
struct FEIterator{FES, IR, G}
    fesp::FES
    _bir::IR
    _geom::G
    _dofs::Vector{Int64}
    _nodes::Vector{Int64}
    _m::Vector{Int64}
    _irs::Vector{IncRel}
    _flds::Vector{FEField}
    
    function FEIterator(fesp::FES) where {FES}
        _bir = baseincrel(fesp.mesh)
        nd = ndofsperelem(fesp.fe)
        _geom = MeshCore.attribute(_bir.right, "geom")
        _dofs = zeros(Int64, nd)
        nn = nfeatofdim(fesp.fe, 0)
        _nodes = zeros(Int64, nn) 
        _m = Int64[]
        _irs = IncRel[]
        _flds = FEField[]
        for m in keys(fesp._irsfields)
            v = fesp._irsfields[m]
            push!(_m, m)
            push!(_irs, v[1])
            push!(_flds, v[2])
        end
        return new{FES, typeof(_bir), typeof(_geom)}(fesp, _bir, _geom, _dofs, _nodes, _m, _irs, _flds)
    end
end

function _storedofs!(d, p, e, ir, fl)
    c = retrieve(ir, e)
    for k in c
        eq = fl.dofnums[k]
        for n in eq
            d[p] = n
            p = p + 1
        end
    end
    return p
end

function update!(it::FEIterator, state)
    #= TODO =#
    copyto!(it._nodes, retrieve(it._bir, state))
    p = 1
    for i in 1:length(it._m)
        m = it._m[i]
        p = _storedofs!(it._dofs, p, state, it._irs[i], it._flds[i])
    end
    return it
end

function Base.iterate(it::FEIterator, state = 1)
    if state > nrelations(it._bir)
        return nothing
    else
        return (update!(it, state), state+1)
    end
end
Base.length(it::FEIterator)  = nrelations(it._bir)

"""
    numberdofs!(self::FEField)

Number the degrees of freedom.

The unknown degrees of freedom in the FE space are numbered consecutively. Then
the datum degrees of freedom are numbered, again consecutively. 

No effort is made to optimize the numbering in any way. 
"""
function numberdofs!(fesp::FES)  where {FES<:FESpace}
    function fieldnumberdofs!(f) 
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

    for m in keys(fesp._irsfields)
        if ndofsperfeat(fesp.fe, m) > 0
            v = fesp._irsfields[m]
            fieldnumberdofs!(v[2]) 
        end 
    end
    return fesp
end

"""
    ndofs(fesp::FES)  where {FES<:FESpace}

Compute the total number of degrees of freedom.
"""
function ndofs(fesp::FES)  where {FES<:FESpace}
    n = 0
    for m in keys(fesp._irsfields)
        if ndofsperfeat(fesp.fe, m) > 0
            v = fesp._irsfields[m]
            n = n + ndofs(v[2]) 
        end 
    end
    return n
end

"""
    nunknowns(fesp::FES)  where {FES<:FESpace}

Compute the total number of unknown degrees of freedom.
"""
function nunknowns(fesp::FES)  where {FES<:FESpace}
    n = 0
    for m in keys(fesp._irsfields)
        if ndofsperfeat(fesp.fe, m) > 0
            v = fesp._irsfields[m]
            n = n + nunknowns(v[2]) 
        end 
    end
    return n
end

"""
    setebc!(fesp::FESpace, mid, eid, comp, val::T) where {T}

Set the EBCs (essential boundary conditions).

- `mid`  = manifold dimension of the entity,
- `eid`  = serial number of the entity (term identifier),
- `comp` = which  degree of freedom in the term,
- `val`  = value of type T

For instance, `mid = 0` means set  the degree of freedom at the vertex `eid`.
"""
function setebc!(fesp::FESpace, mid, eid, comp, val::T) where {T}
    v = fesp._irsfields[mid]
    setebc!(v[2], eid, comp, val)
    return  fesp
end

"""
    gathersysvec(self::FEField{N, T}) where {N, T}

Gather values from the field for the whole system vector.
"""
# function gathersysvec(self::FESpace)
#     nt = nterms(self)
#     ndpt = ndofsperterm(self)
#     vec = fill(zero(eltype(self.dofvals[1])), nt * ndpt)
#     for i in 1:nt
#         en = self.dofnums[i]
#         for j in 1:ndpt
#             vec[en[j]] = self.dofvals[i][j]
#         end
#     end
#     return vec
# end


end

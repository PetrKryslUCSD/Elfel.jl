module FESpaces

using StaticArrays
using MeshCore
using MeshCore: nshapes, indextype, nrelations, retrieve, manifdim
using MeshKeeper: Mesh, baseincrel, increl
using ..FElements: nfeatofdim, ndofsperfeat, ndofsperelem, ndofsperfeat

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

struct FESpace{FET}
    fe::FET
    mesh::Mesh
    _fields::Vector{FEField}

    function FESpace(::Type{T}, fe::FET, mesh) where {T, FET}
        baseir = baseincrel(mesh)
        _fields = _makefields(T, indextype(baseir), fe, mesh)
        return new{FET}(fe, mesh, _fields)
    end
end

function _makefields(::Type{T}, ::Type{IT}, fe, mesh) where {T, IT} 
    fv = [FEField(Val(ndofsperfeat(fe, m)), T, IT, 0) for m in 0:1:3]
    for m in 0:1:manifdim(fe.sd)
        if ndofsperfeat(fe, m) > 0
            ir = increl(mesh, (manifdim(fe.sd), m))
            fv[m+1] = FEField(Val(ndofsperfeat(fe, m)), T, IT, nshapes(ir.right))
        end 
    end
    return fv
end

#= TODO is it more natural to have access to the geometry from the font element space or from the iterator? =#
struct FEIterator{FES, IR, G}
    fesp::FES
    _bir::IR
    _geom::G
    _dofs::Vector{Int64}
    _nodes::Vector{Int64}
    
    function FEIterator(fesp::FES) where {FES}
        _bir = baseincrel(fesp.mesh)
        nd = ndofsperelem(fesp.fe)
        _geom = MeshCore.attribute(_bir.right, "geom")
        _dofs = zeros(Int64, nd)
        nn = nfeatofdim(fesp.fe, 0)
        _nodes = zeros(Int64, nn) 
        return new{FES, typeof(_bir), typeof(_geom)}(fesp, _bir, _geom, _dofs, _nodes)
    end
end

function update!(i::FEIterator, state)
    #= TODO =#
    copyto!(i._nodes, retrieve(i._bir, state))
    return i
end

function Base.iterate(i::FEIterator, state = 1)
    if state > nrelations(i._bir)
        return nothing
    else
        return (update!(i, state), state+1)
    end
end
Base.length(i::FEIterator)  = nrelations(i._bir)

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

    for m in 0:1:manifdim(fesp.fe.sd)
        if ndofsperfeat(fesp.fe, m) > 0
            fieldnumberdofs!(fesp._fields[m+1]) 
        end 
    end
    return fesp
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

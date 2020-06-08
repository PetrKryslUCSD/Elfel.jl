module FESpaces

using StaticArrays
using MeshCore
using MeshCore: nshapes, indextype, nrelations, nentities, retrieve, manifdim, IncRel, VecAttrib
using MeshSteward: Mesh, baseincrel, increl
using ..FElements: nfeatofdim, feathasdof
using ..FEFields: FEField, nterms
import ..FEFields: numberdofs!, ndofs, setebc!, nunknowns, scattersysvec!, gathersysvec!
import ..FElements: ndofsperel

struct FESpace{FET, T}
    mesh::Mesh
    fe::FET
    nfecopies::Int64
    _irsfields::Dict

    function FESpace(::Type{T}, mesh, fe::FET, nfecopies = 1) where {FET, T}
        baseir = baseincrel(mesh)
        _irsfields = _makefields(T, indextype(baseir), mesh, fe, nfecopies)
        return new{FET, T}(mesh, fe, nfecopies, _irsfields)
    end
end

"""
    doftype(fesp::FESpace{FET, T}) where {FET, T}

Provide the type of the values of the degrees of freedom.
"""
doftype(fesp::FESpace{FET, T}) where {FET, T} = T

function _makefields(::Type{T}, ::Type{IT}, mesh, fe, nfecopies) where {T, IT} 
    _irsfields= Dict()
    for m in 0:1:manifdim(fe.sd)
        if feathasdof(fe, m)
            fir = increl(mesh, (manifdim(fe.sd), m))
            fld = FEField(Val(nfecopies), T, IT, nshapes(fir.right))
            _irsfields[m] = (fir, fld)
        end 
    end
    return _irsfields
end

ndofsperel(fesp::FES)  where {FES<:FESpace} = ndofsperel(fesp.fe) * fesp.nfecopies

"""
    numberdofs!(self::FEField)

Number the degrees of freedom.

The unknown degrees of freedom in the FE space are numbered consecutively. Then
the datum degrees of freedom are numbered, again consecutively. 

No effort is made to optimize the numbering in any way. 
"""
function numberdofs!(fesp::FES)  where {FES<:FESpace}
    for m in keys(fesp._irsfields)
        if feathasdof(fesp.fe, m)
            v = fesp._irsfields[m]
            numberdofs!(v[2]) 
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
        if feathasdof(fesp.fe, m) 
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
        if feathasdof(fesp.fe, m)
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
    gathersysvec!(v, fesp::FESpace)

Gather values for the whole system vector.
"""
function gathersysvec!(v, fesp::FESpace)
    for m in keys(fesp._irsfields)
        if feathasdof(fesp.fe, m) 
            gathersysvec!(v, fesp._irsfields[m][2])
        end 
    end
    return v
end

"""
    scattersysvec!(fesp::FESpace, v)

Scatter values from the system vector.
"""
function scattersysvec!(fesp::FESpace, v)
    for m in keys(fesp._irsfields)
        if feathasdof(fesp.fe, m) 
            scattersysvec!(fesp._irsfields[m][2], v)
        end 
    end
    return fesp
end

function makeattribute(fesp::FESpace, name, comp)
    for m in keys(fesp._irsfields)
        if feathasdof(fesp.fe, m)
            ir, fl = fesp._irsfields[m]
            ir.right.attributes[name] = VecAttrib([fl.dofvals[i][comp] for i in 1:nterms(fl)])
        end 
    end
end

end

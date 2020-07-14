module FESpaces

using StaticArrays
using MeshCore
using MeshCore: nshapes, indextype, nrelations, nentities, retrieve, IncRel, VecAttrib
using MeshSteward: Mesh, baseincrel, increl
using ..FElements: nfeatofdim, ndofperfeat, manifdim
import ..FElements: ndofsperel

include("FEFields.jl")

using .FEFields: FEField, nterms
import .FEFields: ndofs, setebc!, scattersysvec!, gathersysvec!
import .FEFields: numberfreedofs!, numberdatadofs!, freedofnums, datadofnums
import .FEFields: highestfreedofnum, highestdatadofnum, dofnums


"""
    FESpace{FET, T}

Type of a finite element space, parameterized with
- `FET`: type of finite element, it is a *scalar* finite element, 
- `T`: type of degree of freedom value (double, complex, ...).
"""
struct FESpace{FET, T}
    mesh::Mesh
    fe::FET
    nfecopies::Int64
    _irsfields::Dict
    _edofmdim::Vector{Int64}
    _edofbfnum::Vector{Int64}
    _edofcompnt::Vector{Int64}

    function FESpace(::Type{T}, mesh, fe::FET, nfecopies = 1) where {FET, T}
        baseir = baseincrel(mesh)
        _irsfields = _makefields(T, indextype(baseir), mesh, fe, nfecopies)
        _edofmdim, _edofbfnum, _edofcompnt = _number_edofs(fe, nfecopies)
        return new{FET, T}(mesh, fe, nfecopies, _irsfields, _edofmdim, _edofbfnum, _edofcompnt)
    end
end

"""
    doftype(fesp::FESpace{FET, T}) where {FET, T}

Provide the type of the values of the degrees of freedom.
"""
doftype(fesp::FESpace{FET, T}) where {FET, T} = T

"""
    edofmdim(fesp::FESpace{FET, T}) where {FET, T}

Access vector of manifold dimensions of entities associated with each degree
of freedom.
"""
edofmdim(fesp::FESpace{FET, T}) where {FET, T} = fesp._edofmdim
"""
    edofbfnum(fesp::FESpace{FET, T}) where {FET, T}

Access vector of numbers of basis functions associated with each degree
of freedom.
"""
edofbfnum(fesp::FESpace{FET, T}) where {FET, T} = fesp._edofbfnum

"""
    edofcompnt(fesp::FESpace{FET, T}) where {FET, T}

Access vector of component number associated with each degree
of freedom.

When the finite element space consists of multiple copies of the scalar finite
element, the component is the serial number of the copy.
"""
edofcompnt(fesp::FESpace{FET, T}) where {FET, T} = fesp._edofcompnt

function _makefields(::Type{T}, ::Type{IT}, mesh, fe, nfecopies) where {T, IT} 
    _irsfields= Dict()
    for m in 0:1:manifdim(fe)
        if ndofperfeat(fe, m) > 0
            fir = increl(mesh, (manifdim(fe), m))
            fld = FEField(Val(nfecopies), T, IT, nshapes(fir.right))
            _irsfields[m] = (fir, fld)
        end 
    end
    return _irsfields
end

function _number_edofs(fe, nfecopies)
    emdim = Int64[]
    bfnum = Int64[]
    compnt = Int64[]
    bfn = 1
    for m in 0:1:3
        for i in 1:nfeatofdim(fe, m) 
            for k in 1:ndofperfeat(fe, m)
                for j in 1:nfecopies
                    push!(emdim, m)
                    push!(bfnum, bfn)
                    push!(compnt, j)
                end
                bfn += 1
            end
        end
    end
    return emdim, bfnum, compnt
end

"""
    ndofsperel(fesp::FES)  where {FES<:FESpace}

Total number of degrees of freedom associated with each finite element.

Essentially a product of the number of the degrees of freedom the scalar
finite element and the number of copies of this element in the space.
"""
ndofsperel(fesp::FES)  where {FES<:FESpace} = ndofsperel(fesp.fe) * fesp.nfecopies

"""
    dofnum(fesp::FES, m, eid)  where {FES<:FESpace}

Provide degree of freedom number for entity `eid` of manifold dimension `m` and
component `comp`.
"""
function dofnum(fesp::FES, m, eid, comp)  where {FES<:FESpace}
    n = zero(typeof(eid))
    if ndofperfeat(fesp.fe, m) > 0
        f = fesp._irsfields[m][2]
        n = dofnums(f, eid)[comp]
    end
    return n
end

"""
    numberfreedofs!(fesp::FES, firstnum = 1)  where {FES<:FESpace}

Number the free degrees of freedom.

The unknown degrees of freedom in the FE space are numbered consecutively. 

No effort is made to optimize the numbering in any way. 
"""
function numberfreedofs!(fesp::FES, firstnum = 1)  where {FES<:FESpace}
    for m in keys(fesp._irsfields)
        if ndofperfeat(fesp.fe, m) > 0
            f = fesp._irsfields[m][2]
            numberfreedofs!(f, firstnum) 
            fnum, lnum, tnum = freedofnums(f)
            firstnum = lnum + 1 
        end 
    end
    return fesp
end

"""
    numberdatadofs!(fesp::FES, firstnum = 1)  where {FES<:FESpace}

Number the data (known) degrees of freedom.

The known degrees of freedom in the FE space are numbered consecutively. 

No effort is made to optimize the numbering in any way. 
"""
function numberdatadofs!(fesp::FES, firstnum = 0)  where {FES<:FESpace}
    firstnum = (firstnum == 0 ? nunknowns(fesp) + 1 : firstnum)
    for m in keys(fesp._irsfields)
        if ndofperfeat(fesp.fe, m) > 0
            f = fesp._irsfields[m][2]
            numberdatadofs!(f, firstnum) 
            fnum, lnum, tnum = datadofnums(f)
            firstnum = lnum + 1 
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
        if ndofperfeat(fesp.fe, m) > 0 
            f = fesp._irsfields[m][2]
            n = n + ndofs(f) 
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
        if ndofperfeat(fesp.fe, m) > 0
            f = fesp._irsfields[m][2]
            fnum, lnum, tnum = freedofnums(f)
            n += tnum
        end 
    end
    return n
end

"""
    highestfreedofnum(fesp::FES)  where {FES<:FESpace}

Compute the highest number of free (unknown) degrees of freedom.
"""
function highestfreedofnum(fesp::FES)  where {FES<:FESpace}
    n = 0
    for m in keys(fesp._irsfields)
        if ndofperfeat(fesp.fe, m) > 0
            f = fesp._irsfields[m][2]
            n = max(n, highestfreedofnum(f))
        end 
    end
    return n
end

"""
    highestfreedofnum(fesp::FES)  where {FES<:FESpace}

Compute the highest number of data (known) degrees of freedom.
"""
function highestdatadofnum(fesp::FES)  where {FES<:FESpace}
    n = 0
    for m in keys(fesp._irsfields)
        if ndofperfeat(fesp.fe, m) > 0
            f = fesp._irsfields[m][2]
            n = max(n, highestdatadofnum(f))
        end 
    end
    return n
end

"""
    numberdofs!(fesp::AbstractVector)

Number the degrees of freedom of a collection of FE spaces.

The unknown (free) degrees of freedom in the FE space are numbered consecutively, and
then the data degrees of freedom (the known values) are numbered. 

No effort is made to optimize the numbering in any way. 
"""
function numberdofs!(fesp::AbstractVector) 
    numberfreedofs!(fesp[1], 1)
    for i in 2:length(fesp)
        numberfreedofs!(fesp[i], highestfreedofnum(fesp[i-1])+1)
    end
    numberdatadofs!(fesp[1], highestfreedofnum(fesp[end])+1)
    for i in 2:length(fesp)
        numberdatadofs!(fesp[i], highestdatadofnum(fesp[i-1])+1)
    end
end

"""
    numberdofs!(fesp::FESpace) 

Number the degrees of freedom of a single FE space.

The unknown (free) degrees of freedom in the FE space are numbered consecutively, and
then the data degrees of freedom (the known values) are numbered. 

No effort is made to optimize the numbering in any way. 
"""
function numberdofs!(fesp::FESpace) 
    numberdofs!([fesp]) 
end

"""
    setebc!(fesp::FESpace, m, eid, comp, val::T) where {T}

Set the EBCs (essential boundary conditions).

- `m`  = manifold dimension of the entity,
- `eid`  = serial number of the entity (term identifier),
- `comp` = which  degree of freedom in the term,
- `val`  = value of type T

For instance, `m = 0` means set  the degree of freedom at the vertex `eid`.
"""
function setebc!(fesp::FESpace, m, eid, comp, val::T) where {T}
    v = fesp._irsfields[m]
    setebc!(v[2], eid, comp, val)
    return  fesp
end

"""
    gathersysvec!(v, fesp::FESpace)

Gather values for the whole system vector.
"""
function gathersysvec!(v, fesp::FESpace)
    for m in keys(fesp._irsfields)
        if ndofperfeat(fesp.fe, m) > 0 
            gathersysvec!(v, fesp._irsfields[m][2])
        end 
    end
    return v
end

"""
    gathersysvec!(v, fesp)

Gather values for the whole system vector from all FE spaces contributing to it.

`fesp` is either a vector or a tuple of FE spaces.
"""
function gathersysvec!(v, fesp::AbstractVector) 
    for i in 1:length(fesp)
        gathersysvec!(v, fesp[i])
    end
    return v
end

"""
    scattersysvec!(fesp::FESpace, v)

Scatter values from the system vector.
"""
function scattersysvec!(fesp::FESpace, v)
    for m in keys(fesp._irsfields)
        if ndofperfeat(fesp.fe, m) > 0 
            scattersysvec!(fesp._irsfields[m][2], v)
        end 
    end
    return fesp
end

"""
    scattersysvec!(fesp, v)

Scatter values for the whole system vector to all FE spaces contributing to it.

`fesp` is either a vector or a tuple of FE spaces.
"""
function scattersysvec!(fesp::AbstractVector, v) 
    for i in 1:length(fesp)
        scattersysvec!(fesp[i], v)
    end
    return fesp
end

"""
    makeattribute(fesp::FESpace, name, comp)

Attach attribute to the right shape collection of all incidence relations. 
"""
function makeattribute(fesp::FESpace, name, comp)
    for m in keys(fesp._irsfields)
        if ndofperfeat(fesp.fe, m) > 0
            ir, fl = fesp._irsfields[m]
            v = fill(SVector{length(comp)}(zeros(length(comp))), nterms(fl))
            for i in 1:nterms(fl)
                v[i] = SVector{length(comp)}(fl.dofvals[i][comp])
            end
            ir.right.attributes[name] = VecAttrib(v)
        end 
    end
end

end

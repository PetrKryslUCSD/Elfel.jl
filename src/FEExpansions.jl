module FEExpansions

using StaticArrays
using MeshCore
using MeshCore: nshapes
using MeshKeeper: Mesh, insert!, increl, basecode
using Elfel.FESpaces: FESpace, multiplicity
using Elfel.Fields: FEField, ndofsperentity, nentities, setebc!, gathersysvec
import Elfel.Fields: numberdofs!, ndofs

"""
    FEExpansion{FET, GT}

Finite element expansion type. 

It links together the mesh, the finite element space, and a field.
The field provides storage of the degree of freedom values and of the numbering.
"""
struct FEExpansion{FET, GT}
    mesh::Mesh
    fesp::FESpace{FET}
    field::FEField
    _geom::GT
end

function FEExpansion(mesh::Mesh, fesp::FESpace{FET}) where {FET}
    irc = basecode(mesh)
    ir = increl(mesh, irc)
    nent = nshapes(ir.right)
    ndpe = multiplicity(fesp)
    f = FEField{ndpe, Float64}(nent)
    _geom = MeshCore.attribute(ir.right, "geom")
    return FEExpansion(mesh, fesp, f, _geom)
end

geometry(fex::FEExpansion) = fex._geom

function numberdofs!(fex::FEExpansion)
    numberdofs!(fex.field)
    return fex
end

ndofs(fex::FEExpansion) = ndofs(fex.field)

end # module

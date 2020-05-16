module FEExpansions

using StaticArrays
using MeshCore
using MeshCore: nshapes
using MeshKeeper: Mesh, insert!, increl, basecode
using Elfel.FESpaces: FESpace, multiplicity
using Elfel.Fields: FEField, ndofsperentity, nentities, setebc!, gathersysvec
import Elfel.Fields: numberdofs!, ndofs

struct FEExpansion{FET, GT}
    mesh::Mesh
    fesp::FESpace{FET}
    field::FEField
    _geomval::GT
end

function FEExpansion(mesh::Mesh, fesp::FESpace{FET}) where {FET}
    irc = basecode(mesh)
    ir = increl(mesh, irc)
    nent = nshapes(ir.right)
    ndpe = multiplicity(fesp)
    f = FEField{ndpe, Float64}(nent)
    _geomval = MeshCore.attribute(ir.right, "geom")
    return FEExpansion(mesh, fesp, f, _geomval)
end

geomval(fex::FEExpansion) = fex._geomval

function numberdofs!(fex::FEExpansion)
    numberdofs!(fex.field)
    return fex
end

ndofs(fex::FEExpansion) = ndofs(fex.field)

end # module

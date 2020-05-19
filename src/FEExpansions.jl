module FEExpansions

using StaticArrays
using MeshCore
using MeshCore: nshapes
using MeshKeeper: Mesh, insert!, increl, basecode
using ..FElements: AbstractFE, refshape
# using Elfel.FESpaces: FESpace, multiplicity
# using Elfel.Fields: FEField, ndofsperentity, nentities, setebc!, gathersysvec
# import Elfel.Fields: numberdofs!, ndofs

"""
    FEExpansion{FET, GT}

Finite element expansion type. 

It links together the mesh, the finite element space, and a field.
The field provides storage of the degree of freedom values and of the numbering.
"""
struct FEExpansion{FET, T} where {FET<:AbstractFE}
    mesh::Mesh
    fe::FET
    _zero::T
end

function FEExpansion(ir::IncRel, fe::FET) where {FET, T}
    baseir = increl(mesh, basecode(mesh))
    nv = nshapes(baseir.right)
    dofnums =  VecAttrib([zero(indextype(baseir)) for i in 1:nv])
    ir.right.attributes["dofnums"] = dofnums
    isdatum =  VecAttrib([false for i in 1:nv])
    ir.right.attributes["isdatum"] = isdatum
    dofvals =  VecAttrib([zero(T) for i in 1:nv])
    ir.right.attributes["dofvals"] = dofvals
    return FEExpansion(mesh, fe, zero(T))
end

"""
    geometry(fex::FEExpansion)

Return the geometry attribute of the base relation of the mesh.
"""
function geometry(fex::FEExpansion)
    ir = increl(mesh, basecode(mesh))
    return MeshCore.attribute(ir.right, "geom")
end

# function numberdofs!(fex::FEExpansion)
#     numberdofs!(fex.field)
#     return fex
# end

# ndofs(fex::FEExpansion) = ndofs(fex.field)

end # module

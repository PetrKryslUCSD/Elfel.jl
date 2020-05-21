module FESpaces

using StaticArrays
using ..RefShapes: RefShapeInterval, RefShapeTriangle
using ..FElements: AbstractFE, refshape
using ..FEMeshes: FEMesh
using ..FEFields: FEField

"""
    FESpace

Type to represent the finite element space.  
"""
struct FESpace
    femesh::FEMesh
    fefield::FEField
end

end # module

module FEMeshes

using StaticArrays
using MeshKeeper: Mesh, baseincrel
using MeshCore: nshapes, indextype, VecAttrib, AbsAttrib, attribute
using ..FElements: nodesperelem, ndofpernode

struct FEMesh
    mesh
    fe
end

geomattr(femesh::FEMesh) = MeshCore.attribute(baseincrel(femesh.mesh).right, "geom")

end

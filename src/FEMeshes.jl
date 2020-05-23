module FEMeshes

using StaticArrays
using MeshKeeper: Mesh, baseincrel
import MeshKeeper: nspacedims
using MeshCore
using MeshCore: nshapes, indextype, VecAttrib, AbsAttrib, attribute, nrelations, retrieve
import ..RefShapes: manifdim
using ..FElements: refshape

struct FEMesh
    mesh
    fe
end

struct FEMeshConn{R}
    ir::R
end

geomattr(femesh::FEMesh) = MeshCore.attribute(baseincrel(femesh.mesh).right, "geom")
nspacedims(femesh::FEMesh) = nspacedims(femesh.mesh)
manifdim(femesh::FEMesh) = manifdim(refshape(femesh.fe))
connectivity(femesh::FEMesh) = FEMeshConn(baseincrel(femesh.mesh))
Base.iterate(e::FEMeshConn, state=1) = begin
    state > nrelations(e.ir) ? 
    nothing : 
    (MeshCore.retrieve(e.ir, state), state+1)
end

end

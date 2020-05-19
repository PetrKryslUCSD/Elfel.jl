module FEMeshes

using StaticArrays
using MeshKeeper: Mesh, baseincrel
using MeshCore: nshapes, indextype, VecAttrib, AbsAttrib
using ..FElements: nodesperelem, ndofpernode

struct FEMesh
    mesh
    fe
    function FEMesh(mesh, fe)
        ir = baseincrel(mesh)
        sc = ir.right
        nv = nshapes(sc)
        npe = nodesperelem(fe)
        ndn = ndofpernode
        @show dofnums =  VecAttrib([zero(typeof(nshapes(sc))) for i in 1:nv])
        @show typeof(dofnums) <: AbsAttrib
        @show typeof(sc.attributes)
        sc.attributes["dofnums"] = dofnums
        isdatum =  VecAttrib([false for i in 1:nv])
        sc.attributes["isdatum"] = isdatum
        dofvals =  VecAttrib([zero(T) for i in 1:nv])
        sc.attributes["dofvals"] = dofvals
        return new(mesh, fe)
    end
end

end

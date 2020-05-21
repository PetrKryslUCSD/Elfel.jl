# module mfesp1
# using Elfel
# using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
# using Elfel.FElements: FE, nodesperelem, refshape
# using Elfel.FElements: bfun, bfundpar, nbasisfuns
# using Elfel.FESpaces: FESpace
# using Test
# function test()
#     e = FEH1_T3(1)
    
#     fesp = FESpace()
#     (e, multiplicity) = fesp.feinfo
#     @test (multiplicity, nbasisfuns(e)) == (1, 3)
        
# end
# end
# using .mfesp1
# mfesp1.test()

module mfes3
using StaticArrays
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, nodesperelem, refshape, FEH1_T3
using Elfel.FElements: bfun, bfundpar
using Elfel.FESpaces: FESpace
using Elfel.FEMeshes: FEMesh, connectivity, iterate, geomattr
using Elfel.FEFields: FEField
using MeshCore
using MeshKeeper: Mesh, load, nspacedims, baseincrel
using Test
function test()
    fe = FEH1_T3(1)
    mesh = load(Mesh(), "qmesh.mesh")
    femesh = FEMesh(mesh, fe)
    fefield = FEField(Float64, femesh::FEMesh) 
    fespace = FESpace(femesh, fefield)
    
    sdim = nspacedims(femesh)
    mdim = manifdim(femesh)

    geom = geomattr(femesh)

    # @show fex.fesp.feinfo[1]
    el = 4
    gradNpar = bfundpar(fespace.femesh.fe, [1/3, 1/3])
    
    conn = connectivity(femesh)
    for c in conn
        J = SMatrix{sdim, mdim}(zeros(sdim, mdim))
        for j in 1:length(c)
            J += geom[c[j]] * gradNpar[j]
        end
    end

    @show c = iterate(conn, 1)[1]
    J = SMatrix{sdim, mdim}(zeros(sdim, mdim))
    for j in 1:length(c)
        J += geom[c[j]] * gradNpar[j]
    end
    @test isapprox(J, [0.5 0.5; 0.0 0.3333333333333333])
end
end
using .mfes3
mfes3.test()

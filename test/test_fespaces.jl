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
using Elfel.FEMeshes: FEMesh
using Elfel.FEFields: FEField
using MeshCore
using MeshKeeper: Mesh, load, nspacedims, baseincrel
using Test
function test()
    fe = FEH1_T3(1)
    mesh = load(Mesh(), "qmesh.mesh")
    femesh = FEMesh(mesh, fe)
    fefield = FEField(Float64, femesh::FEMesh) 
    FEField(Val(ndofpernode(fe)), ::Type{T}, ::Type{IT}, sc::S)
    fespace = FESpace(femesh, fefield)
    fex = FEExpansion(mesh, fesp)
fesp = FESpace((FE{RefShapeTriangle, 3, 1}(), 1))
    
    sdim = nspacedims(fex.mesh)
    ir = baseincrel(fex.mesh)
    mdim = MeshCore.manifdim(ir.left)

    # @test nodesperelem(fex.fesp.feinfo[1]) == 3
    # @test refshape(fex.fesp.feinfo[1]) == RefShapeTriangle
    # @test manifdim(refshape(fex.fesp.feinfo[1])) == 2
    # @test isapprox(bfun(fex.fesp.feinfo[1], [1/3, 1/3]), [0.3333333333333334, 0.3333333333333333, 0.3333333333333333])
    # g = bfundpar(fex.fesp.feinfo[1], [1/3, 1/3])
    # @test isapprox(g[1], [-1.; -1.]')
    # @test isapprox(g[2], [+1.;  0.]')
    # @test isapprox(g[3], [0.; +1.]')

    geom = MeshCore.attribute(ir.right, "geom")

    # @show fex.fesp.feinfo[1]
    el = 4
    gradNpar = bfundpar(fex.fesp.feinfo[1], [1/3, 1/3])
    
    J = SMatrix{sdim, mdim}(zeros(sdim, mdim))
    for j in 1:MeshCore.nentities(ir, 1)
        k = MeshCore.retrieve(ir, 1, j)
        J += geom[k] * gradNpar[j]
    end

    @test isapprox(J, [0.5 0.5; 0.0 0.3333333333333333])

    # e = FE{RefShapeInterval, 2, 1}()
    # @test nodesperelem(e) == 2
    # @test refshape(e) == RefShapeInterval
    # @test manifdim(refshape(e)) == 1
    # @test isapprox(bfun(e, [1/3]), [0.3333333333333334, 2*0.3333333333333333])
    # g = bfundpar(e, [1/3])
    # @test isapprox(g[1], [-1.0/2])
    # @test isapprox(g[2], [+1.0/2])
end
end
using .mfes3
mfes3.test()

module mfes3
using StaticArrays
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, nodesperelem, refshape, FEH1_T3
using Elfel.FElements: bfun, bfundpar
using Elfel.FESpaces: FESpace
using Elfel.FEMeshes: FEMesh
using Elfel.FEFields: FEField
using MeshCore
using MeshKeeper: Mesh, load, nspacedims, baseincrel
using Test
function test()
    fe = FEH1_T3(1)
    mesh = load(Mesh(), "qmesh.mesh")
    femesh = FEMesh(mesh, fe)
    fefield = FEField(Float64, femesh) 
    fespace = FESpace(femesh, fefield)
    
    geom = geomattr(femesh)
    sdim = nspacedims(femesh.mesh)
    mdim = manifdim(refshape(femesh.fe))

    el = 4
    gradNpar = bfundpar(fespace.femesh.fe, [1/3, 1/3])
    
    J = SMatrix{sdim, mdim}(zeros(sdim, mdim))
    for j in 1:MeshCore.nentities(ir, 1)
        k = MeshCore.retrieve(ir, 1, j)
        J += geom[k] * gradNpar[j]
    end

    @test isapprox(J, [0.5 0.5; 0.0 0.3333333333333333])

    # e = FE{RefShapeInterval, 2, 1}()
    # @test nodesperelem(e) == 2
    # @test refshape(e) == RefShapeInterval
    # @test manifdim(refshape(e)) == 1
    # @test isapprox(bfun(e, [1/3]), [0.3333333333333334, 2*0.3333333333333333])
    # g = bfundpar(e, [1/3])
    # @test isapprox(g[1], [-1.0/2])
    # @test isapprox(g[2], [+1.0/2])
end
end
using .mfes3
mfes3.test()
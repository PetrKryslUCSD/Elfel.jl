module mfes1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, nodesperelem, refshape, FEH1_L2, FEH1_T3
using Elfel.FElements: bfun, bfundpar, ndofpernode
using Test
function test()
    e = FE{RefShapeTriangle, 3, 1}()
    e2 = FEH1_T3(1)
    @test e == e2

    @test nodesperelem(e) == 3
    @test refshape(e) == RefShapeTriangle
    @test manifdim(refshape(e)) == 2
    @test isapprox(bfun(e, [1/3, 1/3]), [0.3333333333333334, 0.3333333333333333, 0.3333333333333333])
    g = bfundpar(e, [1/3, 1/3])
    @test isapprox(g[1], [-1.; -1.]')
    @test isapprox(g[2], [+1.;  0.]')
    @test isapprox(g[3], [0.; +1.]')

    e = FE{RefShapeInterval, 2, 1}()
    e2 = FEH1_L2(1)
    @test e == e2
    @test nodesperelem(e) == 2
    @test ndofpernode(e) == 1
    @test refshape(e) == RefShapeInterval
    @test manifdim(refshape(e)) == 1
    @test isapprox(bfun(e, [1/3]), [0.3333333333333334, 2*0.3333333333333333])
    g = bfundpar(e, [1/3])
    @test isapprox(g[1], [-1.0/2])
    @test isapprox(g[2], [+1.0/2])
end
end
using .mfes1
mfes1.test()

# module mfes2
# using Elfel
# using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
# using Elfel.FElements: FE, nodesperelem, refshape
# using Elfel.FElements: bfun, bfundpar
# using Elfel.FESpaces: FESpace
# using MeshKeeper: Mesh, load
# using Elfel.FEExpansions: FEExpansion
# using Test
# function test()
#     fesp = FESpace((FE{RefShapeTriangle, 3}(), 1))
#     mesh = load(Mesh(), "mt3gen3.mesh")
#     fex = FEExpansion(mesh, fesp)

#     e = fex.fesp.feinfo[1]
#     @test nodesperelem(e) == 3
#     @test refshape(e) == RefShapeTriangle
#     @test manifdim(refshape(e)) == 2
#     @test isapprox(bfun(e, [1/3, 1/3]), [0.3333333333333334, 0.3333333333333333, 0.3333333333333333])
#     g = bfundpar(e, [1/3, 1/3])
#     @test isapprox(g[1], [-1.; -1.]')
#     @test isapprox(g[2], [+1.;  0.]')
#     @test isapprox(g[3], [0.; +1.]')

#     e = FE{RefShapeInterval, 2}()
#     @test nodesperelem(e) == 2
#     @test refshape(e) == RefShapeInterval
#     @test manifdim(refshape(e)) == 1
#     @test isapprox(bfun(e, [1/3]), [0.3333333333333334, 2*0.3333333333333333])
#     g = bfundpar(e, [1/3])
#     @test isapprox(g[1], [-1.0/2])
#     @test isapprox(g[2], [+1.0/2])
# end
# end
# using .mfes2
# mfes2.test()
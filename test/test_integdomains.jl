module mfex1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval, npts, param_coords, weights
using Elfel.FElements: FE, nodesperelem, refshape
using Elfel.FElements: bfun, bfundpar, nbasisfuns
using Elfel.Fields: FEField, ndofsperentity, nentities, numberdofs!, setebc!, gathersysvec
using Elfel.FESpaces: FESpace
using Elfel.FEExpansions: FEExpansion, geometry
using MeshKeeper: Mesh, load, baseincrel
using MeshCore: retrieve
using Elfel.IntegDomains: IntegDomain, bfundata, jac, quadrule
using Test
# using BenchmarkTools
function test()
    e = FE{RefShapeTriangle, 3, 1}()

    fesp = FESpace((e, 1))
    mesh = load(Mesh(), "Unit-square-mesh.mesh")
    fex = FEExpansion(mesh, fesp)
    idom = IntegDomain(fex, (kind = :default,))
    @test npts(quadrule(idom)) == 1
    @test isapprox(param_coords(quadrule(idom)), [0.3333333333333333 0.3333333333333333])
    @test isapprox(weights(quadrule(idom)), [0.5])

    Ns, gradNparams = bfundata(idom) 
    @test length(Ns[1]) == 3

    geom = geometry(fex)
    ir = baseincrel(fex.mesh)
    Ns, gradNparams = bfundata(idom) 
    J = jac(geom, retrieve(ir, 1), gradNparams[1])
    # @btime J = jac($geom, retrieve($ir, 1), $gradNparams[1])
    @test isapprox(J, [0.5 0.5; 0.0 0.3333333333333333])
end
end
using .mfex1
mfex1.test()


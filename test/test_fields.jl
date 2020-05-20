module mfld1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, nodesperelem, refshape
using Elfel.FElements: bfun, bfundpar, ndofpernode, FEH1_T3
using MeshKeeper: Mesh, load, increl, baseincrel
using MeshCore: nshapes
using Elfel.FEFields: FEField, nterms, numberdofs!, ndofsperterm, gathersysvec, setebc!
using Test
function test()

    mesh = load(Mesh(), "mt3gen3.mesh")
    fe = FEH1_T3(2)

    fef = FEField(Val(ndofpernode(fe)), Float64, Int64, baseincrel(mesh).right)
    @test nterms(fef) == nshapes(baseincrel(mesh).right)
    @test ndofsperterm(fef) == ndofpernode(fe)
    numberdofs!(fef)
    @test fef.nunknowns == nterms(fef) * ndofsperterm(fef)
    setebc!(fef, 3, 2, 3.0)
    setebc!(fef, 5, 1, 5.0)
    numberdofs!(fef)
    @test fef.nunknowns == nterms(fef) * ndofsperterm(fef) - 2

    v = gathersysvec(fef)
    @test v[nterms(fef) * ndofsperterm(fef) - 1] == 3.0
end
end
using .mfld1
mfld1.test()


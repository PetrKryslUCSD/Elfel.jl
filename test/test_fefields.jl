module mfld1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, ndofsperfeat, refshape
using Elfel.FElements: bfun, bfundpar, FEH1_T3
using MeshKeeper: Mesh, load, increl, baseincrel
using MeshCore: nshapes
using Elfel.FEFields: FEField, nterms, numberdofs!, ndofsperterm, doftype, setebc!, gathersysvec!, ndofs
using Test
function test()

    mesh = load(Mesh(), "mt3gen3.mesh")
    fe = FEH1_T3(2)

    fef = FEField(Val(ndofsperfeat(fe, 0)), Float64, Int64, nshapes(baseincrel(mesh).right))
    @test nterms(fef) == nshapes(baseincrel(mesh).right)
    @test ndofsperterm(fef) == ndofsperfeat(fe, 0)
    numberdofs!(fef)
    @test fef.nunknowns == nterms(fef) * ndofsperterm(fef)
    setebc!(fef, 3, 2, 3.0)
    setebc!(fef, 5, 1, 5.0)
    numberdofs!(fef)
    @test fef.nunknowns == nterms(fef) * ndofsperterm(fef) - 2
    @test doftype(fef) == Float64

    v = fill(zero(doftype(fef)), ndofs(fef))
    gathersysvec!(v, fef)
    @test v[nterms(fef) * ndofsperterm(fef) - 1] == 3.0
end
end
using .mfld1
mfld1.test()


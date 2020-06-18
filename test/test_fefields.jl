module mfld1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape
using Elfel.FElements: bfun, bfungradpar, FEH1_T3
using MeshSteward: Mesh, load, increl, baseincrel
using MeshCore: nshapes
using Elfel.FEFields: FEField, nterms, ndofsperterm, doftype, setebc!, gathersysvec!, ndofs
using Elfel.FEFields: numberfreedofs!, freedofnums, numberdatadofs!, datadofnums
using Test
function test()

    mesh = load(Mesh(), "mt3gen3.mesh")
    fe = FEH1_T3()

    fef = FEField(Val(2), Float64, Int64, nshapes(baseincrel(mesh).right))
    @test nterms(fef) == nshapes(baseincrel(mesh).right)
    @test ndofsperterm(fef) == 2

    # Case of no data degrees of freedom
    numberfreedofs!(fef)
    fnum, lnum, tnum = freedofnums(fef)
    numberdatadofs!(fef, lnum + 1)
    fnum, lnum, tnum = freedofnums(fef)
    @test tnum == nterms(fef) * ndofsperterm(fef)
    fnum, lnum, tnum = datadofnums(fef)
    @test tnum == 0

    # Case of 2 degrees of freedom prescribed
    setebc!(fef, 3, 2, 3.0)
    setebc!(fef, 5, 1, 5.0)
    numberfreedofs!(fef, 1)
    fnum, lnum, tnum = freedofnums(fef)
    numberdatadofs!(fef, lnum + 1)
    fnum, lnum, tnum = freedofnums(fef)
    @test tnum == nterms(fef) * ndofsperterm(fef) - 2
    @test fnum == 1
    @test lnum == nterms(fef) * ndofsperterm(fef) - 2
    fnum, lnum, tnum = datadofnums(fef)
    @test tnum == 2
    @test fnum == nterms(fef) * ndofsperterm(fef) - 2 + 1
    @test lnum == nterms(fef) * ndofsperterm(fef)

    @test doftype(fef) == Float64

    v = fill(zero(doftype(fef)), ndofs(fef))
    gathersysvec!(v, fef)
    @test v[nterms(fef) * ndofsperterm(fef) - 1] == 3.0
end
end
using .mfld1
mfld1.test()


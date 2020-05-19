module mfld1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, nodesperelem, refshape
using Elfel.FElements: bfun, bfundpar, nbasisfuns, FEH1_T3
using MeshKeeper: Mesh, load, increl, baseincrel
using Elfel.FEFields: FEField, nterms
using Test
function test()

    mesh = load(Mesh(), "mt3gen3.mesh")
    fe = FEH1_T3()

    fef = FEField(Float64, baseincrel(mesh).right)
    @test nterms(f) == nshapes(baseincrel(mesh).right)
    
    numberdofs!(f)
    @test f.nunknowns == ne * nd
    setebc!(f, 3, 2, 3.0)
    setebc!(f, 5, 1, 5.0)
    numberdofs!(f)
    @test f.nunknowns == ne * nd - 2

    v = gathersysvec(f)
    @test v[ne * nd - 1] == 3.0
end
end
using .mfld1
mfld1.test()


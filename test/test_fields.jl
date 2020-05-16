module mfld1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, nodesperelem, refshape
using Elfel.FElements: bfun, bfundpar, nbasisfuns
using Elfel.Fields: FEField, ndofsperentity, nentities, numberdofs!, setebc!, gathersysvec
using Test
function test()
    nd = 3
    ne = 5
    f = FEField{nd, Float64}(ne)
    
    @test nentities(f) == ne
    @test ndofsperentity(f) == nd

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


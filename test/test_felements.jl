module mfes1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape, FEH1_L2, FEH1_T3, FEH1_Q4
using Elfel.FElements: bfun, bfundpar, nfeatofdim, ndofsperfeat
using Test
function test()
    e = FEH1_T3(1)
    @test nfeatofdim(e, 0) == 3
    @test refshape(e) == RefShapeTriangle
    @test manifdim(refshape(e)) == 2
    @test ndofsperfeat(e, 0) == 1
    @test isapprox(bfun(e, [1/3, 1/3]), [0.3333333333333334, 0.3333333333333333, 0.3333333333333333])
    g = bfundpar(e, [1/3, 1/3])
    @test isapprox(g[1], [-1.; -1.]')
    @test isapprox(g[2], [+1.;  0.]')
    @test isapprox(g[3], [0.; +1.]')

    e = FEH1_T3(4)
    @test nfeatofdim(e, 0) == 3
    @test refshape(e) == RefShapeTriangle
    @test manifdim(refshape(e)) == 2
    @test ndofsperfeat(e, 0) == 4
    @test isapprox(bfun(e, [1/3, 1/3]), [0.3333333333333334, 0.3333333333333333, 0.3333333333333333])
    g = bfundpar(e, [1/3, 1/3])
    @test isapprox(g[1], [-1.; -1.]')
    @test isapprox(g[2], [+1.;  0.]')
    @test isapprox(g[3], [0.; +1.]')

    e = FEH1_L2(1)
    @test nfeatofdim(e, 0) == 2
    @test ndofsperfeat(e, 0) == 1
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

module mfes2
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval, RefShapeSquare
using Elfel.FElements: FE, refshape, FEH1_L2, FEH1_T3, FEH1_Q4
using Elfel.FElements: bfun, bfundpar, nfeatofdim, ndofsperfeat
using Test
function test()
    e = FEH1_T3(1)
    @test ndofsperfeat(e, 0) == 1
    @test ndofsperfeat(e, 1) == 0
    @test ndofsperfeat(e, 2) == 0
    @test ndofsperfeat(e, 3) == 0
    e = FEH1_T3(3)
    @test ndofsperfeat(e, 0) == 3

end
end
using .mfes2
mfes2.test()
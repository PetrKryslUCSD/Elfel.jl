module mfes1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape, FEH1_L2, FEH1_T3, FEH1_Q4
using Elfel.FElements: bfun, bfungradpar, nfeatofdim, feathasdof
using Test
function test()
    fe = FEH1_T3()
    @test feathasdof(fe, 0) == true
    @test feathasdof(fe, 1) == false
    @test feathasdof(fe, 2) == false
    @test feathasdof(fe, 3) == false

    @test refshape(fe) == RefShapeTriangle
    @test manifdim(refshape(fe)) == 2
    @test nfeatofdim(fe, 0) == 3
    @test nfeatofdim(fe, 1) == 3
    @test nfeatofdim(fe, 2) == 1
    @test nfeatofdim(fe, 3) == 0
    @test isapprox(bfun(fe, [1/3, 1/3]), [0.3333333333333334, 0.3333333333333333, 0.3333333333333333])
    g = bfungradpar(fe, [1/3, 1/3])
    @test isapprox(g[1], [-1.; -1.]')
    @test isapprox(g[2], [+1.;  0.]')
    @test isapprox(g[3], [0.; +1.]')

    e = FEH1_T3()
    @test feathasdof(fe, 0) == true
    @test feathasdof(fe, 1) == false
    @test feathasdof(fe, 2) == false
    @test feathasdof(fe, 3) == false
    @test refshape(fe) == RefShapeTriangle
    @test manifdim(refshape(fe)) == 2
    @test isapprox(bfun(fe, [1/3, 1/3]), [0.3333333333333334, 0.3333333333333333, 0.3333333333333333])
    g = bfungradpar(fe, [1/3, 1/3])
    @test isapprox(g[1], [-1.; -1.]')
    @test isapprox(g[2], [+1.;  0.]')
    @test isapprox(g[3], [0.; +1.]')

    fe = FEH1_L2()
    @test feathasdof(fe, 0) == true
    @test feathasdof(fe, 1) == false
    @test feathasdof(fe, 2) == false
    @test feathasdof(fe, 3) == false
    @test refshape(fe) == RefShapeInterval
    @test manifdim(refshape(fe)) == 1
    @test isapprox(bfun(fe, [1/3]), [0.3333333333333334, 2*0.3333333333333333])
    g = bfungradpar(fe, [1/3])
    @test isapprox(g[1], [-1.0/2])
    @test isapprox(g[2], [+1.0/2])
end
end
using .mfes1
mfes1.test()

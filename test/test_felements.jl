module mfes1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape, FEH1_L2, FEH1_T3, FEH1_Q4
using Elfel.FElements: bfun, bfungradpar, nfeatofdim, ndofperfeat
using Test
function test()
    fe = FEH1_T3()
    @test ndofperfeat(fe, 0) == 1
    @test ndofperfeat(fe, 1) == 0
    @test ndofperfeat(fe, 2) == 0
    @test ndofperfeat(fe, 3) == 0

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
    @test ndofperfeat(fe, 0) == true
    @test ndofperfeat(fe, 1) == 0
    @test ndofperfeat(fe, 2) == 0
    @test ndofperfeat(fe, 3) == 0
    @test refshape(fe) == RefShapeTriangle
    @test manifdim(refshape(fe)) == 2
    @test isapprox(bfun(fe, [1/3, 1/3]), [0.3333333333333334, 0.3333333333333333, 0.3333333333333333])
    g = bfungradpar(fe, [1/3, 1/3])
    @test isapprox(g[1], [-1.; -1.]')
    @test isapprox(g[2], [+1.;  0.]')
    @test isapprox(g[3], [0.; +1.]')

    fe = FEH1_L2()
    @test ndofperfeat(fe, 0) == true
    @test ndofperfeat(fe, 1) == 0
    @test ndofperfeat(fe, 2) == 0
    @test ndofperfeat(fe, 3) == 0
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

module mfes2
using Elfel
using Elfel.RefShapes: manifdim, RefShapeTetrahedron
using Elfel.FElements: FE, refshape, FEH1_T4
using Elfel.FElements: bfun, bfungradpar, nfeatofdim, ndofperfeat
using Test
function test()
    fe = FEH1_T4()
    @test ndofperfeat(fe, 0) == 1
    @test ndofperfeat(fe, 1) == 0
    @test ndofperfeat(fe, 2) == 0
    @test ndofperfeat(fe, 3) == 0

    @test refshape(fe) == RefShapeTetrahedron
    @test manifdim(refshape(fe)) == 3
    @test nfeatofdim(fe, 0) == 4
    @test nfeatofdim(fe, 1) == 6
    @test nfeatofdim(fe, 2) == 4
    @test nfeatofdim(fe, 3) == 1
    @test isapprox(bfun(fe, [1/4, 1/4, 1/4]), [1/4, 1/4, 1/4, 1/4])
    g = bfungradpar(fe, [1/4, 1/4, 1/4])
    @test isapprox(g[1], [-1. -1. -1.])
    @test isapprox(g[2], [+1.  0. 0.])
    @test isapprox(g[3], [0. +1.  0.])
    @test isapprox(g[4], [0. 0. +1. ])
true
end
end
using .mfes2
mfes2.test()

module mfes3
using Elfel
using Elfel.RefShapes: manifdim, RefShapeSquare
using Elfel.FElements: FE, refshape, FEL2_Q4
using Elfel.FElements: bfun, bfungradpar, nfeatofdim, ndofperfeat
using Test
function test()
    fe = FEL2_Q4()
    @test ndofperfeat(fe, 0) == 0
    @test ndofperfeat(fe, 1) == 0
    @test ndofperfeat(fe, 2) == 1
    @test ndofperfeat(fe, 3) == 0

    @test refshape(fe) == RefShapeSquare
    @test manifdim(refshape(fe)) == 2
    @test nfeatofdim(fe, 0) == 4
    @test nfeatofdim(fe, 1) == 4
    @test nfeatofdim(fe, 2) == 1
    @test nfeatofdim(fe, 3) == 0
    @test isapprox(bfun(fe, [1/4, 1/4]), [1.0])
    g = bfungradpar(fe, [1/4, 1/4])
    @test isapprox(g[1], [0. 0. ])
true
end
end
using .mfes3
mfes3.test()

module mfes4
using Elfel
using Elfel.RefShapes: manifdim, RefShapeTetrahedron
using Elfel.FElements: FE, refshape, FEL2_T4
using Elfel.FElements: bfun, bfungradpar, nfeatofdim, ndofperfeat
using Test
function test()
    fe = FEL2_T4()
    @test ndofperfeat(fe, 0) == 0
    @test ndofperfeat(fe, 1) == 0
    @test ndofperfeat(fe, 2) == 0
    @test ndofperfeat(fe, 3) == 1

    @test refshape(fe) == RefShapeTetrahedron
    @test manifdim(refshape(fe)) == 3
    @test nfeatofdim(fe, 0) == 4
    @test nfeatofdim(fe, 1) == 6
    @test nfeatofdim(fe, 2) == 4
    @test nfeatofdim(fe, 3) == 1
    @test isapprox(bfun(fe, [1/4, 1/4, 1/4]), [1.0])
    g = bfungradpar(fe, [1/4, 1/4, 1/4])
    @test isapprox(g[1], [0.0 0.0 0.0])
true
end
end
using .mfes4
mfes4.test()
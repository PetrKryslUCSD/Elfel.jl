module mrs1
using Elfel.RefShapes: RefShapeInterval, manifdim
using Elfel.RefShapes: quadrature, npts, param_coords, weights
using Test
function test()
    @test manifdim(RefShapeInterval) == 1
    qr = quadrature(RefShapeInterval)
    @test npts(qr) == 1
    qr = quadrature(RefShapeInterval, (order = 3,))
    @test npts(qr) == 3
    @test isapprox(param_coords(qr), [-0.774596669241483; 0.0; 0.774596669241483])
    @test isapprox(weights(qr), [0.5555555555555556; 0.888888888888889; 0.5555555555555556])
end
end
using .mrs1
mrs1.test()

module mrs2
using Elfel.RefShapes: RefShapeTriangle, manifdim
using Elfel.RefShapes: quadrature, npts, param_coords, weights
using Test
function test()
    @test manifdim(RefShapeTriangle) == 2
    qr = quadrature(RefShapeTriangle)
    @test npts(qr) == 1
    qr = quadrature(RefShapeTriangle, (npts = 3,))
    @test npts(qr) == 3
    @test isapprox(param_coords(qr), [0.6666666666666666 0.16666666666666666; 0.16666666666666666 0.6666666666666666; 0.16666666666666666 0.16666666666666666])
    @test isapprox(weights(qr), [0.16666666666666666; 0.16666666666666666; 0.16666666666666666])
end
end
using .mrs2
mrs2.test()

module mrs3
using Elfel
using Elfel.RefShapes: RefShapeSquare, manifdim
using Elfel.RefShapes: quadrature, npts, param_coords, weights
using Test
function test()
    @test manifdim(RefShapeSquare) == 2
    qr = quadrature(RefShapeSquare)
    # @show qr
    @test npts(qr) == 1
    qr = quadrature(RefShapeSquare, (order = 2,))
    @test npts(qr) == 4
    # @show param_coords(qr)[3, :]
    @test isapprox(vec(param_coords(qr)[3, :]), [0.577350269189626, -0.577350269189626])
    # @test isapprox(weights(qr), [0.16666666666666666; 0.16666666666666666; 0.16666666666666666])
end
end
using .mrs3
mrs3.test()

module mrs4
using Elfel
using Elfel.RefShapes: RefShapeSquare, manifdim
using Elfel.RefShapes: quadrature, npts, param_coords, weights
using Test
function test()
    @test manifdim(RefShapeSquare) == 2
    qr = quadrature(RefShapeSquare)
    # @show qr
    for o in 1:10
        qr = quadrature(RefShapeSquare, (order = o,))
        @test npts(qr) == o^2
    end
    true
end
end
using .mrs4
mrs4.test()

module mrs5
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim
using Elfel.RefShapes: quadrature, npts, param_coords, weights
using Test
function test()
    @test manifdim(RefShapeTriangle) == 2
   for n in [1, 3, 4, 6, 7, 9, 12, 13]
        qr = quadrature(RefShapeTriangle, (npts = n,))
        @test npts(qr) == n
    end
    true
end
end
using .mrs5
mrs5.test()

module mqpit1
using Elfel
using Elfel.RefShapes: npts, param_coords, weights
using Elfel.FElements: FE, refshape
using Elfel.FElements: bfun, bfungradpar, FEH1_T3
using Elfel.FESpaces: FESpace, ndofs, numberdofs!, setebc!, nunknowns
using MeshKeeper: Mesh, load, baseincrel
using MeshCore: retrieve
using Elfel.QPIterators: QPIterator, bfun, bfungradpar, weight
using Test
# using BenchmarkTools
function test()
    fe = FEH1_T3(1)
    
    qpit = QPIterator(fe, (kind = :default,))
    @test length(qpit) == 1
    # @test isapprox(param_coords(quadrule(idom)), [0.3333333333333333 0.3333333333333333])
    # @test isapprox(weights(quadrule(idom)), [0.5])

    refNs = [[0.3333333333333334, 0.3333333333333333, 0.3333333333333333]]   
    k = 1
    for qp in qpit
        Ns = bfun(qp) 
        @test isapprox(Ns, refNs[k])
        k = k + 1
    end

    qpit = QPIterator(fe, (kind = :default, npts = 3))
    @test length(qpit) == 3
    # @test isapprox(param_coords(quadrule(idom)), [0.3333333333333333 0.3333333333333333])
    # @test isapprox(weights(quadrule(idom)), [0.5])

    refNs = [[0.1666666666666667, 0.6666666666666666, 0.16666666666666666],           
    [0.16666666666666674, 0.16666666666666666, 0.6666666666666666],                   
    [0.6666666666666667, 0.16666666666666666, 0.16666666666666666]]   
    k = 1
    for qp in qpit
        Ns = bfun(qp)  
        @test isapprox(Ns, refNs[k])
        k = k + 1
    end
    k = 1
    for qp in qpit
        Ns = bfun(qp)  
        @test isapprox(Ns, refNs[k])
        k = k + 1
    end

    refgNs = [[-1.0 -1.0], [1.0 0.0], [0.0 1.0]]
    for qp in qpit
        gradNs = bfungradpar(qp) 
        for i in 1:length(gradNs)
            @test isapprox(gradNs[i], refgNs[i])
        end
    end

    for qp in qpit
        w = weight(qp) 
        @test w == 0.16666666666666666
    end
    # Ns, gradNparams = bfundata(idom) 
    # @test length(Ns[1]) == 3

    # geom = geomattr(femesh)
    # conn = connectivity(femesh)
    
    # it = FEIterator(fesp)
    # @show c = iterate(it)
    # Ns, gradNparams = bfundata(idom) 
    # J = jac(geom, c, gradNparams[1])
    # # @btime J = jac($geom, retrieve($ir, 1), $gradNparams[1])
    # @test isapprox(J, [0.5 0.5; 0.0 0.3333333333333333])
end
end
using .mqpit1
mqpit1.test()

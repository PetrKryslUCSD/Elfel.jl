module mfex1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval, npts, param_coords, weights
using Elfel.FElements: FE, nodesperelem, refshape
using Elfel.FElements: bfun, bfundpar, nbasisfuns, FEH1_T3
using Elfel.FEFields: FEField
using Elfel.FEMeshes: FEMesh, connectivity, iterate, geomattr
using MeshKeeper: Mesh, load, baseincrel
using MeshCore: retrieve
using Elfel.IntegDomains: IntegDomain, bfundata, jac, quadrule
using Test
# using BenchmarkTools
function test()
    mesh = load(Mesh(), "Unit-square-mesh.mesh")
    femesh = FEMesh(mesh, FEH1_T3(1))
    fefield = FEField(Float64, femesh) 
    
    idom = IntegDomain(femesh, (kind = :default,))
    @test npts(quadrule(idom)) == 1
    @test isapprox(param_coords(quadrule(idom)), [0.3333333333333333 0.3333333333333333])
    @test isapprox(weights(quadrule(idom)), [0.5])

    Ns, gradNparams = bfundata(idom) 
    @test length(Ns[1]) == 3

    geom = geomattr(femesh)
    conn = connectivity(femesh)
    
    c = iterate(conn, 1)[1]
    Ns, gradNparams = bfundata(idom) 
    J = jac(geom, c, gradNparams[1])
    # @btime J = jac($geom, retrieve($ir, 1), $gradNparams[1])
    @test isapprox(J, [0.5 0.5; 0.0 0.3333333333333333])
end
end
using .mfex1
mfex1.test()

module mfes3
using StaticArrays
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, nodesperelem, refshape, FEH1_T3
using Elfel.FElements: bfun, bfundpar
using Elfel.FEMeshes: FEMesh, connectivity, iterate, geomattr
using Elfel.FEFields: FEField
using MeshCore
using MeshKeeper: Mesh, load, nspacedims, baseincrel
using Test
function test()
    fe = FEH1_T3(1)
    mesh = load(Mesh(), "qmesh.mesh")
    femesh = FEMesh(mesh, fe)
    fefield = FEField(Float64, femesh::FEMesh) 
    
    sdim = nspacedims(femesh)
    mdim = manifdim(femesh)

    geom = geomattr(femesh)

    # @show fex.fesp.feinfo[1]
    el = 4
    gradNpar = bfundpar(femesh.fe, [1/3, 1/3])
    
    conn = connectivity(femesh)
    for c in conn
        J = SMatrix{sdim, mdim}(zeros(sdim, mdim))
        for j in 1:length(c)
            J += geom[c[j]] * gradNpar[j]
        end
    end

    c = iterate(conn, 1)[1]
    J = SMatrix{sdim, mdim}(zeros(sdim, mdim))
    for j in 1:length(c)
        J += geom[c[j]] * gradNpar[j]
    end
    @test isapprox(J, [0.5 0.5; 0.0 0.3333333333333333])
end
end
using .mfes3
mfes3.test()



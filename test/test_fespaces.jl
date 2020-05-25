module mfesp1
using StaticArrays
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape, FEH1_T3
using Elfel.FElements: bfun, bfungradpar
using Elfel.FESpaces: FESpace, ndofs, numberdofs!, setebc!, nunknowns, doftype, nunknowns
using Elfel.FEIterators: FEIterator
using MeshCore
using MeshKeeper: Mesh, load, nspacedims, baseincrel
using Test
function test()
    fe = FEH1_T3(1)
    mesh = load(Mesh(), "qmesh.mesh")
    fesp = FESpace(Float64, fe, mesh)

    @test doftype(fesp) == Float64
    
    refc = [[1, 2, 5],                                                                                                                                                             
    [1, 5, 4],                                                                                                                                                             
    [4, 5, 8],                                                                                                                                                             
    [4, 8, 7],                                                                                                                                                             
    [7, 8, 11],                                                                                                                                                            
    [7, 11, 10],                                                                                                                                                           
    [2, 3, 6],                                                                                                                                                             
    [2, 6, 5],                                                                                                                                                             
    [5, 6, 9],                                                                                                                                                             
    [5, 9, 8],                                                                                                                                                             
    [8, 9, 12],                                                                                                                                                            
    [8, 12, 11], ]   
    it = FEIterator(fesp)
    k = 1
    for c in it
        @test isapprox(c._nodes,  refc[k])
        k = k + 1
    end
    @test length(it) == 12

    # @show summary(mesh)

    numberdofs!(fesp)
    @test ndofs(fesp) == 12

    for i in [1, 4, 7, 10]
        setebc!(fesp, 0, i, 1, 0.0)
    end
    numberdofs!(fesp)
    @test ndofs(fesp) == 12
    @test nunknowns(fesp) == 8

    it = FEIterator(fesp)
    @time for el in it
         el._dofs
    end
    # sdim = nspacedims(femesh)
    # mdim = manifdim(femesh)

    # geom = geomattr(femesh)  

    # el = 4
    # gradNpar = bfungradpar(fespace.fe, [1/3, 1/3])
    
    # conn = connectivity(femesh)
    # for c in conn
    #     J = SMatrix{sdim, mdim}(zeros(sdim, mdim))
    #     for j in 1:length(c)
    #         J += geom[c[j]] * gradNpar[j]
    #     end
    # end

    # c = iterate(conn, 1)[1]
    # J = SMatrix{sdim, mdim}(zeros(sdim, mdim))
    # for j in 1:length(c)
    #     J += geom[c[j]] * gradNpar[j]
    # end
    # @test isapprox(J, [0.5 0.5; 0.0 0.3333333333333333])
end
end
using .mfesp1
mfesp1.test()

# module mfes3
# using StaticArrays
# using Elfel
# using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
# using Elfel.FElements: FE, refshape, FEH1_T3
# using Elfel.FElements: bfun, bfungradpar
# using Elfel.FESpaces: FESpace
# using MeshCore
# using MeshKeeper: Mesh, load, nspacedims, baseincrel
# using Test
# function test()
#     fe = FEH1_T3(1)
#     mesh = load(Mesh(), "qmesh.mesh")
#     fesp = FESpace(mesh, fe)
    
#     sdim = nspacedims(femesh)
#     mdim = manifdim(femesh)

#     geom = geomattr(femesh)

#     el = 4
#     gradNpar = bfungradpar(fespace.fe, [1/3, 1/3])
    
#     conn = connectivity(femesh)
#     for c in conn
#         J = SMatrix{sdim, mdim}(zeros(sdim, mdim))
#         for j in 1:length(c)
#             J += geom[c[j]] * gradNpar[j]
#         end
#     end

#     c = iterate(conn, 1)[1]
#     J = SMatrix{sdim, mdim}(zeros(sdim, mdim))
#     for j in 1:length(c)
#         J += geom[c[j]] * gradNpar[j]
#     end
#     @test isapprox(J, [0.5 0.5; 0.0 0.3333333333333333])
# end
# end
# using .mfes3
# mfes3.test()

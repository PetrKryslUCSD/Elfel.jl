module maggfeit1
using StaticArrays
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape, FEH1_T3
using Elfel.FElements: bfun, bfungradpar
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype
using Elfel.FEIterators: FEIterator, elnodes
using Elfel.AggregateFEIterators: AggregateFEIterator, iterator
using Elfel.Utilities: @_check

using MeshCore
using MeshSteward: Mesh, load, nspacedims, baseincrel
using Test
function test()
    mesh = load(Mesh(), "qmesh.mesh")
    fesp1 = FESpace(Float64, mesh, FEH1_T3())
    fesp2 = FESpace(Float64, mesh, FEH1_T3())
    fesp3 = FESpace(Float64, mesh, FEH1_T3())

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
    it1 = FEIterator(fesp1)
    it2 = FEIterator(fesp2)
    it3 = FEIterator(fesp3)
    it = AggregateFEIterator([it1, it2, it3])
    @test length(it) == 12

    k = 1
    for i in it
        @test elnodes(iterator(i, 1)) == elnodes(iterator(i, 2)) == elnodes(iterator(i, 3))
        k = k + 1
    end
    
    
end
end
using .maggfeit1
maggfeit1.test()


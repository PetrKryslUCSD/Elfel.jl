module mfeit1
using StaticArrays
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape, FEH1_T3
using Elfel.FElements: bfun, bfungradpar
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype
using Elfel.FESpaces: numberfreedofs!, numberdatadofs!
using Elfel.FEIterators: FEIterator

using MeshCore
using MeshSteward: Mesh, load, nspacedims, baseincrel
using Test
function test()
    fe = FEH1_T3()
    mesh = load(Mesh(), "qmesh.mesh")
    fesp = FESpace(Float64, mesh, fe)

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

    numberfreedofs!(fesp)
    numberdatadofs!(fesp)
    @test ndofs(fesp) == 12
    @test nunknowns(fesp) == 12

    for i in [1, 4, 7, 10]
        setebc!(fesp, 0, i, 1, 0.0)
    end
    numberfreedofs!(fesp)
    numberdatadofs!(fesp)
    @test ndofs(fesp) == 12
    @test nunknowns(fesp) == 8

    # @show fesp._irsfields

    refd =  [[9], [1], [2], [10], [3], [4], [11], [5], [6], [12], [7], [8]]
    it = FEIterator(fesp)
    for el in it
         # @show el._dofs
         # @show refd[el._nodes]
         @test isapprox(el._dofs, [refd[n][1] for n in el._nodes] )
    end
end
end
using .mfeit1
mfeit1.test()

module mfeit2
using StaticArrays
using LinearAlgebra
using MeshCore
using MeshSteward: Mesh, load, nspacedims, baseincrel
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape, FEH1_T3
using Elfel.FElements: bfun, bfungradpar
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype
using Elfel.FESpaces: numberfreedofs!, numberdatadofs!
using Elfel.FEIterators: FEIterator, eldofs
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Elfel.LocalAssemblers: LocalMatrixAssembler, LocalVectorAssembler, init!
using Test
A = [0.6744582963441466 0.2853043149927861 0.27460710155821255; 
0.3781923479141225 0.2838873430062512 0.6316949656630075; 
0.19369805365903336 0.8926164783344779 0.07006905962860177]
function test()
    fe = FEH1_T3()
    mesh = load(Mesh(), "qmesh.mesh")
    fesp = FESpace(Float64, mesh, fe)

    for i in [1, 4, 7, 10]
        setebc!(fesp, 0, i, 1, 0.0)
    end
    numberfreedofs!(fesp)
    numberdatadofs!(fesp)

    ass = SysmatAssemblerSparse(0.0)
    start!(ass, 12, 12)
    ke = LocalMatrixAssembler(size(A, 1), size(A, 2), 0.0)
    it = FEIterator(fesp)
    for el in it
        init!(ke, eldofs(el), eldofs(el))
        for j in 1:size(A, 2), i in 1:size(A, 1)
            ke[i, j] += A[i, j]
        end
        assemble!(ass, ke)
    end
    S = finish!(ass)

    D = fill(0.0, 12, 12)
    it = FEIterator(fesp)
    for el in it
        for j in 1:size(A, 2), i in 1:size(A, 1)
            D[el._dofs[i], el._dofs[j]] += A[i, j]
        end
    end
    @test isapprox(D - S, 0 * I)
end
end
using .mfeit2
mfeit2.test()

module mfeit3
using StaticArrays
using LinearAlgebra
using MeshCore
using MeshSteward: Mesh, load, nspacedims, baseincrel
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape, FEH1_T3
using Elfel.FElements: bfun, bfungradpar
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype
using Elfel.FESpaces: numberfreedofs!, numberdatadofs!
using Elfel.FEIterators: FEIterator, eldofs, eldofs, eldofentmdims, eldofcomps
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Test
A = [0.6744582963441466 0.2853043149927861 0.27460710155821255; 
0.3781923479141225 0.2838873430062512 0.6316949656630075; 
0.19369805365903336 0.8926164783344779 0.07006905962860177]
function test()
    fe = FEH1_T3()
    mesh = load(Mesh(), "qmesh.mesh")
    fesp = FESpace(Float64, mesh, fe, 3)

    for i in [1, 4, 7, 10]
        setebc!(fesp, 0, i, 1, 0.0)
    end
    numberfreedofs!(fesp)
        numberdatadofs!(fesp)

    it = FEIterator(fesp)
    @test em = eldofentmdims(it) == [0, 0, 0, 0, 0, 0, 0, 0, 0]  
    @test eldofcomps(it) == [1, 2, 3, 1, 2, 3, 1, 2, 3] 
    ref = [ [33, 1, 2, 3, 4, 5, 11, 12, 13],                              
     [33, 1, 2, 11, 12, 13, 34, 9, 10],                            
     [34, 9, 10, 11, 12, 13, 19, 20, 21],                          
     [34, 9, 10, 19, 20, 21, 35, 17, 18],                          
     [35, 17, 18, 19, 20, 21, 27, 28, 29],                         
     [35, 17, 18, 27, 28, 29, 36, 25, 26],                         
     [3, 4, 5, 6, 7, 8, 14, 15, 16],                               
     [3, 4, 5, 14, 15, 16, 11, 12, 13],                            
     [11, 12, 13, 14, 15, 16, 22, 23, 24],                         
     [11, 12, 13, 22, 23, 24, 19, 20, 21],                         
     [19, 20, 21, 22, 23, 24, 30, 31, 32],                         
     [19, 20, 21, 30, 31, 32, 27, 28, 29]  ]
    i = 1
    for el in it
        @test isapprox(eldofs(it), ref[i])
        i = i + 1
    end

end
end
using .mfeit3
mfeit3.test()

module maggfeit1
using StaticArrays
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape, FEH1_T3
using Elfel.FElements: bfun, bfungradpar
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype
using Elfel.FEIterators: FEIterator, elnodes
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
    it = zip(it1, it2, it3)
    @test length(it) == 12

    k = 1
    for i in it
        @test elnodes(i[1]) == elnodes(i[2]) == elnodes(i[3])
        k = k + 1
    end
    
    
end
end
using .maggfeit1
maggfeit1.test()
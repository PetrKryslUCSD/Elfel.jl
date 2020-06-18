module mfesp1
using StaticArrays
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape, FEH1_T3
using Elfel.FElements: bfun, bfungradpar
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype, nunknowns
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
    @test ndofs(fesp) == 12
    @test nunknowns(fesp) == 12

    for i in [1, 4, 7, 10]
        setebc!(fesp, 0, i, 1, 0.0)
    end
    numberfreedofs!(fesp)
    @test ndofs(fesp) == 12
    @test nunknowns(fesp) == 8

    it = FEIterator(fesp)
    for el in it
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

module mfesp2
using StaticArrays
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape, FEH1_T3, shapedesc
using Elfel.FElements: FEH1_T3_BUBBLE
using Elfel.FElements: bfun, bfungradpar
using Elfel.FElements: nfeatofdim, ndofperfeat, ndofsperel
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype, nunknowns
using Elfel.FESpaces: numberfreedofs!, numberdatadofs!
using Elfel.FESpaces: edofmdim, edofbfnum, edofcompnt
using Elfel.FEIterators: FEIterator
using MeshCore: identty
using MeshSteward: Mesh, load, nspacedims, baseincrel, insert!
using Test
function test()
    fe = FEH1_T3()
    mesh = load(Mesh(), "qmesh.mesh")
    fesp = FESpace(Float64, mesh, fe, 3)

    @test doftype(fesp) == Float64
    
    @test shapedesc(fesp.fe).name == "T3"
    @test refshape(fesp.fe) == Elfel.RefShapes.RefShapeTriangle

    @test nfeatofdim(fesp.fe, 0) == 3
    @test ndofperfeat(fesp.fe, 0) == 1

    @test ndofsperel(fesp.fe) == 3

    @test isapprox(bfun(fesp.fe, [1/3, 1/3]), [0.3333333333333334, 0.3333333333333333, 0.3333333333333333])

    insert!(mesh, identty(baseincrel(mesh)))
    fesp = FESpace(Float64, mesh, FEH1_T3_BUBBLE(), 1)

    @test nfeatofdim(fesp.fe, 0) == 3
    @test nfeatofdim(fesp.fe, 1) == 3
    @test nfeatofdim(fesp.fe, 2) == 1
    @test nfeatofdim(fesp.fe, 3) == 0
    @test ndofperfeat(fesp.fe, 0) == 1
    @test ndofperfeat(fesp.fe, 1) == 0
    @test ndofperfeat(fesp.fe, 2) == 1
    @test ndofperfeat(fesp.fe, 3) == 0

    @test ndofsperel(fesp.fe) == 4

    @test isapprox(bfun(fesp.fe, [1/3, 1/3]), [0.3333333333333334, 0.3333333333333333, 0.3333333333333333, 0.03703703703703704])

    fesp = FESpace(Float64, mesh, FEH1_T3_BUBBLE(), 2)

    @test nfeatofdim(fesp.fe, 0) == 3
    @test nfeatofdim(fesp.fe, 1) == 3
    @test nfeatofdim(fesp.fe, 2) == 1
    @test nfeatofdim(fesp.fe, 3) == 0
    @test ndofperfeat(fesp.fe, 0) == 1
    @test ndofperfeat(fesp.fe, 1) == 0
    @test ndofperfeat(fesp.fe, 2) == 1
    @test ndofperfeat(fesp.fe, 3) == 0

    @test ndofsperel(fesp.fe) == 4

    @test isapprox(bfun(fesp.fe, [1/3, 1/3]), [0.3333333333333334, 0.3333333333333333, 0.3333333333333333, 0.03703703703703704])

    emdim = Int64[]
    bfnum = Int64[]
    compnt = Int64[]
    bfn = 1
    for m in 0:1:3
        if ndofperfeat(fesp.fe, m) > 0
            for i in 1:nfeatofdim(fesp.fe, m) 
                for k in 1:ndofperfeat(fesp.fe, m)
                    for j in 1:fesp.nfecopies
                        push!(emdim, m)
                        push!(bfnum, bfn)
                        push!(compnt, j)
                    end
                    bfn += 1
                end
            end
        end
    end
    @test edofmdim(fesp) == emdim
    @test edofbfnum(fesp) == bfnum
    @test edofcompnt(fesp) == compnt

    true
end
end
using .mfesp2
mfesp2.test()


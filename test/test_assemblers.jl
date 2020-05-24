module mvass1
using Elfel.Assemblers: SysvecAssembler, start!, finish!, assemble!
using LinearAlgebra
using Test
function test()
    N = 30
    v = fill(0.0, N)
    a = SysvecAssembler()
    start!(a,  N)
    vec, dofnums = rand(3), [1, 7, 2]
    for i in 1:length(vec)
        assemble!(a, dofnums[i], vec[i]) 
    end
    for (i, p) in zip(dofnums, vec)
        v[i] += p
    end
    vec, dofnums = rand(7), [29, 15, 1, 7, 3, 6, 2]
    for i in 1:length(vec)
        assemble!(a, dofnums[i], vec[i]) 
    end
    for (i, p) in zip(dofnums, vec)
        v[i] += p
    end
    w = finish!(a)
    @test norm(v-w) / norm(w) <= 1.0e-9
end
end
using .mvass1
mvass1.test()

module mass1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, refshape
using Elfel.FElements: bfun, bfundpar
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Test
function test()
    a = SysmatAssemblerSparse(0.0)                                                        
    start!(a, 7, 7)  
    m = [0.24406   0.599773    0.833404  0.0420141                                             
    0.786024  0.00206713  0.995379  0.780298                                              
    0.845816  0.198459    0.355149  0.224996]     
    gi = [1 7 5]             
    gj = [5 2 1 4]       
    for j in 1:size(m, 2), i in 1:size(m, 1)
          assemble!(a, gi[i], gj[j], m[i, j])       
    end  

    m = [0.146618  0.53471   0.614342    0.737833                                              
    0.479719  0.41354   0.00760941  0.836455                                              
    0.254868  0.476189  0.460794    0.00919633                                            
    0.159064  0.261821  0.317078    0.77646                                               
    0.643538  0.429817  0.59788     0.958909]                                   
    gi =  [2 3 1 7 5]
    gj = [6 7 3 4]   
    for j in 1:size(m, 2), i in 1:size(m, 1)
          assemble!(a, gi[i], gj[j], m[i, j])       
    end                               
    A = finish!(a) 
    @test abs.(maximum([0.833404 0.0 0.460794 0.05121043 0.24406 0.254868 0.476189; 
        0.0 0.0 0.614342 0.737833 0.0 0.146618 0.53471; 0.0 0.0 0.00760941 0.836455 0.0 0.479719 0.41354; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.355149 0.0 0.59788 1.183905 0.845816 0.643538 0.429817; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.995379 0.0 0.0 0.780298 0.786024 0.0 0.0]  - A)) < 1.0e-5
    #     a = SysmatAssemblerSparse(0.0)
    #     startassembly!(a, 5, 5, 3, 7, 7)
    #     m = [0.24406   0.599773    0.833404  0.0420141
    #     0.786024  0.00206713  0.995379  0.780298
    #     0.845816  0.198459    0.355149  0.224996]
    #     assemble!(a, m, [1 7 5], [5 0 1 4])
    #     m = [0.146618  0.53471   0.614342    0.737833
    #  0.479719  0.41354   0.00760941  0.836455
    #  0.254868  0.476189  0.460794    0.00919633
    #  0.159064  0.261821  0.317078    0.77646
    #  0.643538  0.429817  0.59788     0.958909]
    #     assemble!(a, m, [2 3 1 0 5], [6 7 3 4])
    #     A = makematrix!(a)
    #     # @show Matrix(A)
    #     @test abs.(maximum([0.833404 0.0 0.460794 0.05121043 0.24406 0.254868 0.476189; 
    # 0.0 0.0 0.614342 0.737833 0.0 0.146618 0.53471; 0.0 0.0 0.00760941 0.836455 0.0 0.479719 0.41354; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.355149 0.0 0.59788 1.183905 0.845816 0.643538 0.429817; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.995379 0.0 0.0 0.780298 0.786024 0.0 0.0]  - A)) < 1.0e-5
        
end
end
using .mass1
mass1.test()

module massa1
using Elfel.Assemblers: LocalMatrixAssembler, initialize!, assemble!
using Test
function test()
	A = [0.6744582963441466 0.2853043149927861 0.27460710155821255; 
	0.3781923479141225 0.2838873430062512 0.6316949656630075; 
	0.19369805365903336 0.8926164783344779 0.07006905962860177]
	la = LocalMatrixAssembler(size(A, 1), size(A, 2), 0.0)
	initialize!(la, i -> i, [1, 2, 3])
	# @show rs, cs
	@test isapprox(la.row, [1, 2, 3, 1, 2, 3, 1, 2, 3])
	@test isapprox(la.col, [1, 1, 1, 2, 2, 2, 3, 3, 3])
	for j in 1:size(A, 2), i in 1:size(A, 1)
		assemble!(la, i, j, A[i, j])
	end
	@test isapprox(la.M, A)
	# @show A, vs
end
end
using .massa1
massa1.test()


module massa2
using Elfel.Assemblers: LocalMatrixAssembler, initialize!, assemble!
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Test
function test()
    A = [0.6744582963441466 0.2853043149927861 0.27460710155821255; 
    0.3781923479141225 0.2838873430062512 0.6316949656630075; 
    0.19369805365903336 0.8926164783344779 0.07006905962860177]
    la = LocalMatrixAssembler(size(A, 1), size(A, 2), 0.0)
    initialize!(la, i -> i, [1, 2, 3])
    # @show rs, cs
    @test isapprox(la.row, [1, 2, 3, 1, 2, 3, 1, 2, 3])
    @test isapprox(la.col, [1, 1, 1, 2, 2, 2, 3, 3, 3])
    for j in 1:size(A, 2), i in 1:size(A, 1)
        assemble!(la, i, j, A[i, j])
    end
    ass = SysmatAssemblerSparse(0.0)
    start!(ass, 3, 3)
    assemble!(ass, la)
    As =  finish!(ass)
    matched = 0
    for j in 1:size(A, 2), i in 1:size(A, 1)
       isapprox(A[i, j], As[i, j]) && (matched = matched + 1)
    end
    @test matched == prod(size(A))
    # @show A, vs
end
end
using .massa2
massa2.test()
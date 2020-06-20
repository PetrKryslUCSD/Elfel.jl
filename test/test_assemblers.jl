module mmass0
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using LinearAlgebra
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
    B = [ 0.833404  0.599773    0.460794    0.0512104  0.24406   0.254868  0.476189
    0         0          0.614342    0.737833    0        0.146618  0.53471
    0         0          0.00760941  0.836455    0        0.479719  0.41354
    0         0           0           0          0         0         0
    0.355149  0.198459    0.59788     1.1839     0.845816  0.643538  0.429817
    0         0           0           0          0         0         0
    0.995379  0.00206713  0.317078    1.55676    0.786024  0.159064  0.261821  ]
    @test norm(A - B) / norm(B) < 1.0e-5
    true
end
end
using .mmass0
mmass0.test()

module mmass0a
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Elfel.LocalAssemblers: LocalMatrixAssembler, LocalVectorAssembler, init!
using LinearAlgebra
using Test
function test()
    a = SysmatAssemblerSparse(0.0)                                                        
    start!(a, 7, 7)  
    m = [0.24406   0.599773    0.833404  0.0420141                                             
    0.786024  0.00206713  0.995379  0.780298                                              
    0.845816  0.198459    0.355149  0.224996]     
    gi = [1 7 5]             
    gj = [5 2 1 4]       
    la = LocalMatrixAssembler(3, 4, 0.0)
    init!(la, gi, gj)
    for j in 1:size(m, 2), i in 1:size(m, 1)
        la[i, j] += m[i, j]
    end  
    assemble!(a, la)
    m = [0.146618  0.53471   0.614342    0.737833                                              
    0.479719  0.41354   0.00760941  0.836455                                              
    0.254868  0.476189  0.460794    0.00919633                                            
    0.159064  0.261821  0.317078    0.77646                                               
    0.643538  0.429817  0.59788     0.958909]                                   
    gi =  [2 3 1 7 5]
    gj = [6 7 3 4]   
    la = LocalMatrixAssembler(5, 4, 0.0)
    init!(la, gi, gj)
    for j in 1:size(m, 2), i in 1:size(m, 1)
        la[i, j] += m[i, j]   
    end    
    assemble!(a, la)                  
    A = finish!(a) 
    B = [ 0.833404  0.599773    0.460794    0.0512104  0.24406   0.254868  0.476189
    0         0          0.614342    0.737833    0        0.146618  0.53471
    0         0          0.00760941  0.836455    0        0.479719  0.41354
    0         0           0           0          0         0         0
    0.355149  0.198459    0.59788     1.1839     0.845816  0.643538  0.429817
    0         0           0           0          0         0         0
    0.995379  0.00206713  0.317078    1.55676    0.786024  0.159064  0.261821  ]
    @test norm(A - B) / norm(B) < 1.0e-5
    true
end
end
using .mmass0a
mmass0a.test()

module mmass0t
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Elfel.LocalAssemblers: LocalMatrixAssembler, LocalVectorAssembler, init!
using LinearAlgebra
using Test
function test()
    a = SysmatAssemblerSparse(0.0)                                                        
    start!(a, 7, 7)  
    m = transpose([0.24406   0.599773    0.833404  0.0420141                                             
    0.786024  0.00206713  0.995379  0.780298                                              
    0.845816  0.198459    0.355149  0.224996])     
    gj = [1 7 5]             
    gi = [5 2 1 4]       
    la = LocalMatrixAssembler(4, 3, 0.0)
    @test size(la) == (4, 3)
    init!(la, gi, gj)
    for j in 1:size(m, 2), i in 1:size(m, 1)
        la[i, j] += m[i, j]
    end  
    assemble!(a, transpose(la))
    m = transpose([0.146618  0.53471   0.614342    0.737833                                              
    0.479719  0.41354   0.00760941  0.836455                                              
    0.254868  0.476189  0.460794    0.00919633                                            
    0.159064  0.261821  0.317078    0.77646                                               
    0.643538  0.429817  0.59788     0.958909]   )                                
    gj =  [2 3 1 7 5]
    gi = [6 7 3 4]   
    la = LocalMatrixAssembler(4, 5, 0.0)
    init!(la, gi, gj)
    for j in 1:size(m, 2), i in 1:size(m, 1)
        la[i, j] += m[i, j]   
    end    
    assemble!(a, transpose(la))
    A = finish!(a) 
    B = [ 0.833404  0.599773    0.460794    0.0512104  0.24406   0.254868  0.476189
    0         0          0.614342    0.737833    0        0.146618  0.53471
    0         0          0.00760941  0.836455    0        0.479719  0.41354
    0         0           0           0          0         0         0
    0.355149  0.198459    0.59788     1.1839     0.845816  0.643538  0.429817
    0         0           0           0          0         0         0
    0.995379  0.00206713  0.317078    1.55676    0.786024  0.159064  0.261821  ]
    @test norm(A - B) / norm(B) < 1.0e-5
    la = LocalVectorAssembler(4, 0.0)
    @test size(la)  == (4, )
    true
end
end
using .mmass0t
mmass0t.test()

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
using Elfel.FElements: bfun, bfungradpar
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

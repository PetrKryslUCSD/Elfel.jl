using Test

@time @testset "Assemblers" begin include("test_assemblers.jl") end
@time @testset "Reference shapes" begin include("test_refshapes.jl") end
@time @testset "Finite elements" begin include("test_felements.jl") end
@time @testset "Fields" begin include("test_fefields.jl") end
@time @testset "Spaces" begin include("test_fespaces.jl") end
@time @testset "Quadrature iterators" begin include("test_qpiterators.jl") end
@time @testset "FE iterators" begin include("test_feiterators.jl") end
@time @testset "Heat" begin include("test_heat.jl") end
@time @testset "Elasticity" begin include("test_elasticity.jl") end
@time @testset "Stokes" begin include("test_stokes.jl") end



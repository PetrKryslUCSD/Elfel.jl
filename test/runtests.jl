using Test

@time @testset "Fields" begin include("test_fefields.jl") end
@time @testset "Spaces" begin include("test_fespaces.jl") end
@time @testset "Iterators" begin include("test_feiterators.jl") end
@time @testset "Integration domains" begin include("test_integdomains.jl") end
@time @testset "Reference shapes" begin include("test_refshapes.jl") end
@time @testset "Finite elements" begin include("test_felements.jl") end
@time @testset "Assemblers" begin include("test_assemblers.jl") end


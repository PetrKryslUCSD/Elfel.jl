using Test

@time @testset "Integration domains" begin include("test_integdomains.jl") end
@time @testset "Assemblers" begin include("test_assemblers.jl") end
@time @testset "Fields" begin include("test_fields.jl") end
@time @testset "Reference shapes" begin include("test_refshapes.jl") end
@time @testset "Finite elements" begin include("test_felements.jl") end
@time @testset "FE Spaces" begin include("test_fespaces.jl") end
@time @testset "Expansions" begin include("test_feexpansions.jl") end


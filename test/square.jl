@testset "CIS square example with $optimizer_constructor" for optimizer_constructor in blp_factories
    include(dirname(dirname(pathof(HybridSystems))), "examples", "square.jl")
    horizontal_jump_example(Polyhedra.DefaultLibrary{Float64}(), false)
    p = getis(hs, optimizer_constructor)
    @test p[1].c ≈ [0, 0] atol=1e-5
    @test p[1].Q ≈ [1 0; 0 1] rtol=1e-5
end

@testset "CIS square example with $factory" for factory in blp_factories
    include(dirname(dirname(pathof(HybridSystems))), "examples", "square.jl")
    horizontal_jump_example(Polyhedra.DefaultLibrary{Float64}(), false)
    p = getis(hs, factory)
    @test p[1].c ≈ [0, 0] atol=1e-5
    @test p[1].Q ≈ [1 0; 0 1] rtol=1e-5
end

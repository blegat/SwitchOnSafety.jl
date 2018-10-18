@testset "CIS horizontal jump example with $factory" for factory in blp_factories
    include(dirname(dirname(pathof(HybridSystems))), "examples", "square.jl")
    hs = square_example(Polyhedra.DefaultLibrary{Float64}())
    p = getis(hs, factory)
    @test p[1].c ≈ [0, 0] atol=1e-5
    @test p[1].Q ≈ [1 0; 0 1] rtol=1e-5
    @test p[2].c ≈ [3.14286, 0]
    @test p[2].Q ≈ [1.36111 0
                    0       1.0]
end

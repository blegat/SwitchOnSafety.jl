@testset "CIS horizontal jump example with $solver" for solver in blp_solvers
    include(Pkg.dir("HybridSystems", "examples", "horizontal_jump.jl"));
    p = getis(hs, solver)
    @test p[1].P ≈ [-1 0 0; 0 1 0; 0 0 1]
    @test p[2].P ≈ [168 -70 0; -70 28 0; 0 0 7.0]
end

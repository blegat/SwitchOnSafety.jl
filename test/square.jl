@testset "CIS square example with $solver" for solver in blp_solvers
    include(Pkg.dir("HybridSystems", "examples", "square.jl"))
    p = getis(hs, solver)
    @test p[1].c ≈ [0, 0] atol=1e-5
    @test p[1].Q ≈ [1 0; 0 1] rtol=1e-5
end

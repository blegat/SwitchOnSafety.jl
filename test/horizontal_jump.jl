@testset "CIS horizontal jump example with $solver" for solver in blp_solvers
    include(Pkg.dir("HybridSystems", "examples", "horizontal_jump.jl"));
    p = getis(hs, solver)
    @test p[1].P ≈ [-1 0 0; 0 1 0; 0 0 1] rtol=1e-5
    display(p[2].P)
    @test p[2].P ≈ [128.011  -43.9945  0.0
                    -43.9945  14.0028  0.0
                      0.0      0.0    10.2857] rtol=1e-5
end

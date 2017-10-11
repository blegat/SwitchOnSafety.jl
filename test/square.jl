@testset "CIS square example with $solver" for solver in blp_solvers
    include(Pkg.dir("HybridSystems", "examples", "square.jl"))
    p = getis(hs, solver, [[0., 0]])
    @test p[1].P ≈ [-1 0 0; 0 1 0; 0 0 1]
end

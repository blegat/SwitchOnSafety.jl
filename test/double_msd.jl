using LinearAlgebra

@testset "Quick bounds update smp" begin
    Sigma = [Matrix(1.0I, 4, 4) .+ 1e-4]
    s = discreteswitchedsystem(Sigma)

    qlb, qub = quickb(s, 1)

    @test qlb > 1.

    @test qub > 1.

    @test hassmp(s)

    @test getlb(s) > 1.

    @test getub(s) > 1.
end

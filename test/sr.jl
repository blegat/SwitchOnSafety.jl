@testset "Spectral radius" begin
    @test_throws ArgumentError ρ(zeros(0, 0))
end

using StaticArrays
@testset "Spectral radius" begin
    @test_throws ArgumentError ρ(zeros(0, 0))
    # There is a specialized function for StaticArrays
    @test ρ(@SMatrix [1 2; 3 4]) ≈ ρ([1 2; 3 4])
    @test ρ(-@SMatrix [1 2; 3 4]) ≈ ρ(-[1 2; 3 4])
    @test ρ(@SMatrix [1 1; 0 1]) ≈ ρ([1 1; 0 1])
end

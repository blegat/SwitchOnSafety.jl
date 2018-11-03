using StaticArrays
@testset "Spectral radius" begin
    @test_throws ArgumentError ρ(zeros(0, 0))
    # There is a specialized function for StaticArrays
    @test ρ(@SMatrix [1 2; 3 4]) ≈ ρ([1 2; 3 4])
    @test ρ(-@SMatrix [1 2; 3 4]) ≈ ρ(-[1 2; 3 4])
    @test ρ(@SMatrix [1 1; 0 1]) ≈ ρ([1 1; 0 1])
    # Matrix A1 from [Section 4, PEDJ16]
    A = @SMatrix [0.94  0.56; -0.35  0.73]
    @test ρ(A) ≈ ρ(Matrix(A))
    @test ρ(-A) ≈ ρ(Matrix(-A))
end

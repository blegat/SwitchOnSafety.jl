using Combinatorics

@testset "Veronese Lift" begin
    @test_throws DimensionMismatch permanent([1 2], collect(Combinatorics.permutations(1:2)))
    @test_throws DimensionMismatch SwitchOnSafety.veroneselift_explicit([1 2], 3)
    @test_throws ArgumentError SwitchOnSafety.veroneselift_explicit([1 2; 3 4], 0)
    @test_throws ArgumentError veroneselift([1 2; 3 4], 0)

    @test [1.0 2 * √3 4 * √3 8.0] ≈ @inferred veroneselift([1 2], 3)

    s3 = sqrt(3)
    expected = [11^3 s3*11^2*12 s3*11*12^2 12^3;
    s3*11^2*21 11*(11*22+2*21*12) 12*(2*11*22+21*12) s3*12^2*22;
    s3*11*21^2 21*(2*11*22+21*12) 22*(11*22+2*21*12) sqrt(3)*12*22^2;
    21^3 s3*21^2*22 s3*21*22^2 22^3]
    @test expected ≈ @inferred veroneselift([11 12; 21 22], 3)
    @test expected ≈ @inferred SwitchOnSafety.veroneselift_explicit([11 12; 21 22], 3)
end

# Example 2.8 of
# P. Parrilo and A. Jadbabaie
# Approximation of the joint spectral radius using sum of squares
# Linear Algebra and its Applications, Elsevier, 2008, 428, 2385-2402
# The JSR was conjectured to be 8.9149 = √ρ(A_1 * A_3)

using LinearAlgebra

@testset "[PJ08] Example 5.4" begin
    # values with log-accuracy 4e-7 taken from the examples/PJ08e54.ipynb
    expected_lb = [5.635326998733677, 6.777596245727604, 7.423337986847027]
    expected_ub = [9.760675006197351, 8.91982041593713,  8.914964296278484]
    A1 = [ 0  1  7  4;
           1  6 -2 -3;
          -1 -1 -2 -6;
           3  0  9  1]
    A2 = [-3  3  0 -2;
          -2  1  4  9;
           4 -3  1  1;
           1 -5 -1 -2]
    A3 = [ 1  4  5 10;
           0  5  1 -4;
           0 -1  4  6;
          -1  5  0  1]
    s = discreteswitchedsystem([A1, A2, A3])
    smp = periodicswitching(s, [1, 3])

    @testset "Gripenberg" begin
        for p in [1, 2, Inf]
            psw, ub = gripenberg(s, max_length = 1, matrix_norm = A -> opnorm(A, p), verbose=0)
            @test psw == periodicswitching(s, [2])
            @test ub ≈ opnorm(A3, p)
        end
        psw, ub = gripenberg(s, verbose=0)
        @test psw == smp
        @test ub ≈ psw.growthrate + 1e-2 # 1e-2 is the default δ
    end

    for factory in sdp_factories
        iscsdp(factory) && continue # Segfault on Travis: https://travis-ci.org/blegat/SwitchOnSafety.jl/jobs/588075484
        sosdata(s).lb = 0
        for d in 1:3
            tol = 1e-5
            lb, ub = soslyapb(s, d, factory=factory, tol=tol)
            @test abs(log(expected_lb[d]) - log(lb)) <= tol
            @test abs(log(expected_ub[d]) - log(ub)) <= tol
            seq = sosbuildsequence(s, d, p_0=:Primal)
            psw = findsmp(seq)
            @test psw !== nothing
            @test psw == smp
        end
    end
end

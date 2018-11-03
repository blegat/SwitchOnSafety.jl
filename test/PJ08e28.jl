# Example 2.8 of
# P. Parrilo and A. Jadbabaie
# Approximation of the joint spectral radius using sum of squares
# Linear Algebra and its Applications, Elsevier, 2008, 428, 2385-2402
# Inspired from the construction of:
# Ando, T. and Shih, M.-h.
# Simultaneous Contractibility.
# SIAM Journal on Matrix Analysis & Applications, 1998, 19, 487
# The JSR is 1

@testset "[PJ08] Example 2.8" begin
    expected_lb = [1, 1/2^(1/4)]
    expected_ub = [√2, 1]
    for factory in sdp_factories
        s = discreteswitchedsystem([[1 0; 1 0], [0 1; 0 -1]])
        tol = ismosek(factory) ? 2e-4 : 1e-3
        for d in 1:2
            # We set lb=0 for d=2, otherwise there might be no infeasibile problem and hence no extraction done
            sosdata(s).lb = 0
            lb, ub = soslyapb(s, d, factory=factory, tol=tol)
            @test abs(log(expected_lb[d]) - log(lb)) <= tol
            @test abs(log(expected_ub[d]) - log(ub)) <= tol

            if d == 1
                @test !hassmp(s)
                @test_throws ErrorException getsmp(s)
            else
                @test getlb(s) ≈ 1
                @test hassmp(s)
                @test getsmp(s).growthrate == 1
            end
        end

        for d in 1:2
            seq = sosbuildsequence(s, d, p_0=:Primal)
            psw = findsmp(seq)
            @test psw !== nothing
            @test psw.growthrate == 1
        end
    end
end

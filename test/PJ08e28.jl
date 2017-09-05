# Example 2.8 of
# P. Parrilo and A. Jadbabaie
# Approximation of the joint spectral radius using sum of squares
# Linear Algebra and its Applications, Elsevier, 2008, 428, 2385-2402
# Inspired from the construction of:
# Ando, T. and Shih, M.-h.
# Simultaneous Contractibility.
# SIAM Journal on Matrix Analysis & Applications, 1998, 19, 487
# The JSR is √2

@testset "[PJ08] Example 2.8" begin
    expected_lb = [1, 1/2^(1/4)]
    expected_ub = [√2, 1]
    for solver in sdp_solvers
        s = DiscreteSwitchedSystem([[1 0; 1 0], [0 1; 0 -1]])
        println("  > With solver $(typeof(solver))")
        tol = ismosek(solver) ? 2e-4 : 1e-3
        for d in 1:2
            lb, ub = soslyapb(s, d, solver=solver, tol=tol)
            @test abs(log(expected_lb[d]) - log(lb)) <= tol
            @test abs(log(expected_ub[d]) - log(ub)) <= tol

            if d == 1
                @test isnull(s.smp)
            else
                @test s.lb ≈ 1
                @test !isnull(s.smp)
                @test get(s.smp).growthrate == 1
            end
        end

        for d in 1:2
            seq = sosbuildsequence(s, d, p_0=:Primal)
            psw = findsmp(seq)
            @test !isnull(psw)
            @test get(psw).growthrate == 1
        end
    end
end

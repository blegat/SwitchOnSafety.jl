# Example 2.8 of
# P. Parrilo and A. Jadbabaie
# Approximation of the joint spectral radius using sum of squares
# Linear Algebra and its Applications, Elsevier, 2008, 428, 2385-2402
# The JSR was conjectured to be 8.9149 = √ρ(A_1 * A_3)

@testset "[PJ08] Example 5.4" begin
    expected_lb = [5.635095912315392, 6.777216272557359, 8.92 / 3^(1/6)]
    expected_ub = [9.761463937881695, 8.922585166517319, 8.92]
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
    s = DiscreteSwitchedSystem([A1, A2, A3])
    smp = DiscretePeriodicSwitching(s, [1, 3])
    for solver in sdp_solvers
        s.lb = 0
        iscsdp(solver) && continue
        println("  > With solver $(typeof(solver))")
        for d in 1:3
            tol = ismosek(solver) ? (d <= 2 ? 4e-4 : 3e-2) : 1e-3
            lb, ub = soslyapb(s, d, solver=solver, tol=tol)
            @test abs(log(expected_lb[d]) - log(lb)) <= tol
            @test abs(log(expected_ub[d]) - log(ub)) <= tol
            seq = sosbuildsequence(s, d, p_0=:Primal)
            psw = findsmp(seq)
            @test !isnull(psw)
            @test get(psw) == smp
        end
    end
end

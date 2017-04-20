# Example 6.1 of
# R. Jungers, A. Cicone and N. Guglielmi,
# "Lifted polytope methods for computing the joint spectral radius."
# SIAM Journal on Matrix Analysis and Applications, SIAM, 2014, 35, 391-410.
# The JSR is 1.7779

@testset "[JCG14] Example 6.1" begin
    expected_lb = [1.2626, 1.4949]
    expected_ub = [1.7857, 1.7779]
    s = DiscreteSwitchedSystem([[ 0 -1  1  1;
                                  1  0  0  0;
                                  0 -1  0  0;
                                  1 -1 -1  0],
                                [ 0 -1  1  0;
                                 -1 -1  1  1;
                                 -1  0  0  0
                                 -1 -1  0 -1]])
    for solver in sdp_solvers
        s.lb = 0
        println("  > With solver $(typeof(solver))")
        tol = 1e-4
        for d in 1:2
            lb, ub = soslyapb(s, d, solver=solver, tol=tol)
            @test isapprox(lb, expected_lb[d], rtol=tol)
            @test isapprox(ub, expected_ub[d], rtol=tol)
            seq = sosbuildsequence(s, 1, p_0=:Primal)
            psw = findsmp(seq)
            @test !isnull(psw)
            @test get(psw).period == [2]
            @test isapprox(get(psw).growthrate, 1.7779191220330814)
        end
    end
end

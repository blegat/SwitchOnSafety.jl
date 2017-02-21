# Example 2.1 of
# A. Ahmadi, and P. Parrilo
# "Joint spectral radius of rank one matrices and the maximum cycle mean problem."
# CDC, 731-733, 2012.
# The JSR is 1

@testset "[AP12] Example 2.1" begin
    s = DiscreteSwitchedSystem([[0 1; 0 0], [0 0; 1 0]])
    for solver in sdp_solvers
        println("  > With solver $(typeof(solver))")
        tol = ismosek(solver) ? 1e-5 : 1e-4
        lb, ub = soslyapb(s, 1, solver=solver, tol=tol)
        @test isapprox(lb, 1 / √2, rtol=tol)
        @test isapprox(ub, 1, rtol=tol)
        seq = sosbuildsequence(s, 1, p_0=:Primal)
        psw = findsmp(seq)
        @test !isnull(psw)
        @test get(psw).growthrate == 1
    end
end

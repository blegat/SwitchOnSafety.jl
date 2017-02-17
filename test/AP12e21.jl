# Example 2.1 of
# A. Ahmadi, and P. Parrilo
# "Joint spectral radius of rank one matrices and the maximum cycle mean problem."
# CDC, 731-733, 2012.
# The JSR is 1

@testset "[AP12] Example 2.1" begin
    s = DiscreteSwitchedSystem([0 1; 0 0], [0 0; 1 0])
    for solver in sdp_solvers
        tol = 1e-5
        println("  > With solver $(typeof(solver))")
        lb, ub = soslyapb(s, 1, solver=solver, tol=tol)
        @show lb
        @show ub
    end
end

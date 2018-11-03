# Example 6.1 of
# R. Jungers, A. Cicone and N. Guglielmi,
# "Lifted polytope methods for computing the joint spectral radius."
# SIAM Journal on Matrix Analysis and Applications, SIAM, 2014, 35, 391-410.
# The JSR is 1.7779

@testset "[JCG14] Example 6.1" begin
    expected_lb = [1.2626, 1.4949]
    expected_ub = [1.7857, 1.7779]
    s = discreteswitchedsystem([[ 0 -1  1  1;
                                  1  0  0  0;
                                  0 -1  0  0;
                                  1 -1 -1  0],
                                [ 0 -1  1  0;
                                 -1 -1  1  1;
                                 -1  0  0  0
                                 -1 -1  0 -1]])
    smp = periodicswitching(s, [2])
    for factory in sdp_factories
        sosdata(s).lb = 0
        tol = 1e-4
        for d in 1:2
            lb, ub = soslyapb(s, d, factory=factory, tol=tol)
            @test isapprox(lb, expected_lb[d], rtol=tol)
            @test isapprox(ub, expected_ub[d], rtol=tol)
            seq = sosbuildsequence(s, d, p_0=:Primal)
            psw = findsmp(seq)
            @test psw !== nothing
            @test psw == smp
        end
    end
end

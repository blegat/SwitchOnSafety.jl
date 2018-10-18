# Example 2.1 of
# A. Ahmadi, and P. Parrilo
# "Joint spectral radius of rank one matrices and the maximum cycle mean problem."
# CDC, 731-733, 2012.
# The JSR is 1

@testset "[AP12] Example 2.1" begin
    for factory in sdp_factories
        tol = 1e-4
        for d in 1:3
            # Get a fresh system to discard lyapunovs, smp and sets lb to 0
            s = discreteswitchedsystem([[0 1; 0 0], [0 0; 1 0]])
            smp = periodicswitching(s, [1, 2])
            @test smp.growthrate == 1
            @test periodicswitching(s, [1, 1]) != periodicswitching(s, [2, 2])

            lb, ub = soslyapb(s, d, factory=factory, tol=tol)
            @test lb ≈ 1 / 2^(1/(2*d)) rtol=tol
            @test ub ≈ 1 rtol=tol

            @test getlb(s) == 1.
            @test hassmp(s)
            @test getsmp(s) == smp

            cyc = sosextractcycle(s, d, tol=tol)
            @test cyc !== nothing
            @test cyc == smp

            seq = sosbuildsequence(s, d, p_0=:Primal)
            psw = findsmp(seq)
            @test psw !== nothing
            @test psw.growthrate == 1
        end
    end
end

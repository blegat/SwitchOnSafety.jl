# Section 4 of
# M. Philippe, R. Essick, G. E. Dullerud and R. M. Jungers.
# "Stability of discrete-time switching systems with constrained switching sequences."
# Automatica, 72:242-250, 2016
# The JSR is 0.9748171979372074

using LightGraphs

const expected_ub = [1.0035568400762171,
                     0.9863261446338583,
                     0.9769402921375085,
                     0.9749706065652288,
                     0.9749659545861287,
                     0.9749138871770363]
const ρAs = 2.618033988749896
const ratio = [2, ρAs, ρAs, ρAs, ρAs, ρAs].^(1./(2:2:12))
const expected_lb = expected_ub ./ ratio

@testset "[PEDJ] Section 4" begin
    include(Pkg.dir("HybridSystems", "examples", "PEDJ16s4.jl"))
    sbp = periodicswitching(hs, Edge.([3 => 3]))
    @test sbp.growthrate == 0.9392550239418472
    snp = periodicswitching(hs, Edge.([3 => 1, 1 => 3, 3 => 1, 1 => 3, 3 => 3, 3 => 3, 3 => 3, 3 => 3]))
    @test snp.growthrate == 0.9728940109399586
    smp = periodicswitching(hs, Edge.([3 => 1, 1 => 3, 3 => 1, 1 => 2, 2 => 3, 3 => 3, 3 => 3, 3 => 3]))
    @test smp.growthrate == 0.9748171979372074
    for solver in sdp_solvers
        sosdata(hs).lb = 0
        println("  > With solver $(typeof(solver))")
        for d in 1:6
            tol = ismosek(solver) ? 4e-4 : 1e-3
            sosdata(hs).lb = 0
            lb, ub = soslyapb(hs, d, solver=solver, tol=tol)
            @test log(lb) ≈ log(expected_lb[d]) atol=tol
            @test log(ub) ≈ log(expected_ub[d]) atol=tol
#            sosextractcycle(s, d)
            for l in 1:2
                for v_0 in 1:4
                    seq = sosbuildsequence(hs, d, p_0=:Primal, v_0=v_0)
                    psw = findsmp(seq)
                    @test !isnull(psw)
                    if d <= 3
                        if 2 <= d && v_0 == 3
                            @test get(psw) == sbp
                        else
                            @test get(psw) == snp
                        end
                    else
                        @test get(psw) == smp
                    end
                end
            end
        end
    end
end

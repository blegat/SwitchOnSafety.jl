# Section 4 of
# M. Philippe, R. Essick, G. E. Dullerud and R. M. Jungers.
# "Stability of discrete-time switching systems with constrained switching sequences."
# Automatica, 72:242-250, 2016
# The JSR is 0.9748171979372074

using LightGraphs

const expected_ub = [1.0035568400762171 0.9820958128927598
                     0.9863261446338583 0.9753566329765574
                     0.9769402921375085 0.9751185378641121
                     0.9749706065652288 0.9751185378641121
                     0.9749659545861287 0.9749374081223462
                     0.9749138871770363 0.9748393637263851]
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
    hsm = mdependentlift(hs, 2)
    msbp = periodicswitching(hsm, Edge.([5 => 5]))
    @test msbp.growthrate == 0.9392550239418472
    msnp = periodicswitching(hsm, Edge.([5 => 5, 5 => 3, 3 => 8, 8 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 3, 3 => 8, 8 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 3, 3 => 8, 8 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 3, 3 => 8, 8 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 5]))
    @test msnp.growthrate == 0.9653214971459174
    msmp = periodicswitching(hsm, Edge.([5 => 3, 3 => 8, 8 => 3, 3 => 7, 7 => 2, 2 => 5, 5 => 5, 5 => 5]))
    @test msmp.growthrate == 0.9748171979372074
    for solver in sdp_solvers
        println("  > With solver $(typeof(solver))")
        for s in (hs, hsm)
            m = s === hs ? 1 : 2
            for d in 1:(7-m)
                tol = ismosek(solver) ? 6e-4 : 1e-3
                lb, ub = soslyapb(s, d, solver=solver, tol=tol)
                @test log(lb) ≈ log(expected_lb[d, m]) atol=tol
                @test log(ub) ≈ log(expected_ub[d, m]) atol=tol
            end
            for d in 1:(7-m)
                for l in 1:2
                    for v_0 in states(s)
                        seq = sosbuildsequence(s, d, p_0=:Primal, v_0=v_0, niter=100)
                        psw = findsmp(seq)
                        @test !isnull(psw)
                        if s === hsm
                            if d == 1
                                if v_0 == 5
                                    @test get(psw) == msbp
                                else
                                    @test get(psw) == msnp
                                end
                            elseif d == 2 && v_0 in [3, 5, 8]
                                @test get(psw) == msbp
                            else
                                @test get(psw) == msmp
                            end
                        elseif d <= 3
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
end

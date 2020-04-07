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
const ratio = [2, ρAs, ρAs, ρAs, ρAs, ρAs] .^ (1 ./ (2:2:12))
const expected_lb = expected_ub ./ ratio

@testset "[PEDJ] Section 4" begin
    include(joinpath(dirname(dirname(pathof(HybridSystems))), "examples", "PEDJ16s4.jl"))

    @testset "p-radius with $algo" for algo in (VeroneseLift(), KroneckerLift())
        pρlb = [0.87907456299, 0.89372555286, 0.90488811163]
        pρub = [1.42237252156, 1.13683646451, 1.06232506675]
        for d in 1:3
            lb, ub = pradiusb(hs, 2d, algo)
            @test isapprox(lb, pρlb[d])
            @test isapprox(ub, pρub[d])
        end
    end

    function t(pair::Pair)
        edge = HybridSystems.edge_object(hs.automaton, pair.first, pair.second)
        id = first(keys(hs.automaton.Σ[edge]))
        return HybridSystems.LightTransition(edge, id)
    end
    sbp = periodicswitching(hs, t.([3 => 3]))
    @test sbp.growthrate == 0.9392550239418471
    sap = periodicswitching(hs, t.([3 => 1, 1 => 3, 3 => 1, 1 => 3, 3 => 1, 1 => 3, 3 => 3, 3 => 3, 3 => 3,
                                    3 => 1, 1 => 3, 3 => 1, 1 => 3, 3 => 1, 1 => 3, 3 => 3, 3 => 3, 3 => 3, 3 => 3]))
    @test sap.growthrate == 0.96076337693942
    snp = periodicswitching(hs, t.([3 => 1, 1 => 3, 3 => 1, 1 => 3, 3 => 3, 3 => 3, 3 => 3, 3 => 3]))
    @test snp.growthrate == 0.9728940109399586
    smp = periodicswitching(hs, t.([3 => 1, 1 => 3, 3 => 1, 1 => 2, 2 => 3, 3 => 3, 3 => 3, 3 => 3]))
    @test smp.growthrate == 0.9748171979372074

    hsm = mdependentlift(hs, 2)
    function t(pair::Pair)
        edge = Edge(pair)
        id = first(keys(hsm.automaton.Σ[edge]))
        return HybridSystems.LightTransition(edge, id)
    end
    msbp = periodicswitching(hsm, t.([5 => 5]))
    @test msbp.growthrate == 0.9392550239418471
    msnp = periodicswitching(hsm, t.([5 => 5, 5 => 3, 3 => 8, 8 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 3, 3 => 8, 8 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 3, 3 => 8, 8 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 3, 3 => 8, 8 => 5, 5 => 5, 5 => 5, 5 => 5, 5 => 5]))
    @test msnp.growthrate == 0.9653214971459174
    msmp = periodicswitching(hsm, t.([5 => 3, 3 => 8, 8 => 3, 3 => 7, 7 => 2, 2 => 5, 5 => 5, 5 => 5]))
    @test msmp.growthrate == 0.9748171979372074

    @testset "Polytope with $optimizer_constructor" for optimizer_constructor in soc_optimizer_constructors
        psw, isemp, sets = invariant_polytopes(hs, optimizer_constructor, sbp, tol=1e-4, max_length=8, verbose=0)
        @test psw == smp
        @test isemp
    end

    @testset "Polyset with $optimizer_constructor" for optimizer_constructor in sdp_factories
        @testset "With $(s == hs ? "original system" : "2-dependent lift")" for s in (hs, hsm)
            m = s === hs ? 1 : 2
            @testset "SOS d=$d" for d in 1:(7-m)
                tol = ismosek(optimizer_constructor) ? 6e-4 : 1e-3
                lb, ub = soslyapb(s, d, optimizer_constructor=optimizer_constructor, tol=tol)
                @test log(lb) ≈ log(expected_lb[d, m]) atol=tol
                @test log(ub) ≈ log(expected_ub[d, m]) atol=tol
            end
            @testset "EXT d=$d" for d in 1:(7-m)
                @testset "l=$l" for l in 1:2
                    @testset "v_0=$v_0" for v_0 in states(s)
                        seq = sosbuildsequence(s, d, p_0=:Primal, v_0=v_0, niter=100)
                        psw = findsmp(seq)
                        @test psw !== nothing
                        if s === hsm
                            if d == 1
                                if v_0 == 5
                                    @test psw == msbp
                                else
                                    @test psw == msnp
                                end
                            else
                                # I initially got only msmp but Mosek now also gets msbp
                                @test psw == msmp || psw == msbp
                            end
                        elseif d <= 3
                            if 2 <= d && v_0 == 3
                                @test psw == sbp
                            elseif d == 1 && v_0 in [3, 4]
                                # See https://travis-ci.org/blegat/SwitchOnSafety.jl/jobs/655778174#L568-L665
                                @test psw == sap || psw == snp
                            else
                                @test psw == snp
                            end
                        else
                            @test psw == smp
                        end
                    end
                end
            end
        end
    end
end

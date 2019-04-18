# Example 5.4 of
# A. Ahmadi, R. Jungers, P. Parrilo and M. Roozbehani,
# "Joint spectral radius and path-complete graph Lyapunov functions."
# SIAM J. CONTROL OPTIM 52(1), 687-717, 2014.
# The JSR is 3.917384715148

@testset "[AJPR14] Example 5.4" begin
    s = discreteswitchedsystem([[-1 -1; -4 0], [3 3; -2 1]])
    @test sprint(show, s) == "Hybrid System with automaton OneStateAutomaton(2)"
    qub = 4.31959610746
    @test quicklb(s) ≈ 3
    @test quickub(s) ≈ qub
    qb = quickb(s)
    @test typeof(qb) == NTuple{2,Float64}
    @test qb[1] ≈ 3
    @test qb[2] ≈ qub
    # pnorm=Inf: [17.2409,7.42377,5.65853,4.92192,4.52597,4.28099,4.11275,3.99139,3.8991,3.827,3.76891,3.72121,3.6813,3.64744,3.61834,3.59307,3.57092,3.55135,3.53393,3.51832,3.50426,3.49152,3.47994]
    # pnorm=1:   [15.1245,7.35044,5.47402,4.84722,4.45847,4.23242,4.0709,3.95644,3.8687,3.80015,3.74485,3.69944,3.66141,3.62914,3.6014,3.5773,3.55616,3.53748,3.52085,3.50595,3.49253,3.48037,3.4693,3.45919]
    # pnorm=2:   [14.1535,6.81871,5.32651,4.7096,4.37067,4.15723,4.01117,3.90488,3.82407,3.7606,3.70945,3.66735,3.6321,3.60215,3.57639,3.554,3.53437,3.517,3.50154,3.48768,3.47519]

    @test pradius(s, 2) ≈ 3.234535151244 rtol=1e-9
    @test pradius(s, 2, BruteForce(), pnorm=2) ≈ 3.632097274822649

    pρlb = [3.23453515151, 3.42756015606, 3.54374928704]
    pρub = [4.57432347886, 4.07607892469, 3.97772408342]
    @testset "p-radius with $algo" for algo in (VeroneseLift(), KroneckerLift())
        for d in 1:3
            lb, ub = pradiusb(s, 2d, algo)
            @test isapprox(lb, pρlb[d])
            @test isapprox(ub, pρub[d])
        end
    end

    tmp = periodicswitching(s, [1, 2])
    tmp = periodicswitching(tmp.s, tmp.period, round(tmp.growthrate, digits=4)) # Avoid difference of floating point rounding
    @test sprint(show, tmp) == "PSW(3.9174, [1, 2])" #"Periodic switching of growth rate 3.9174 and modes: [1, 2] for the transitions [1, 2]"
    smp = periodicswitching(s, [1, 2])
    #@test smp != tmp
    @test smp != periodicswitching(s, [2, 1, 2, 1])
    @test smp != periodicswitching(s, [1, 2, 1, 2])
    @test smp != periodicswitching(s, [1, 1])
    #@test smp != periodicswitching(discreteswitchedsystem([[-1 -1; -4 0],[3 3; -2 1]]), [1, 2])
    @testset "JSR with $factory" for factory in sdp_factories
        sosdata(s).lb = 0
        tol = ismosek(factory) ? 1e-5 : 5e-4
        lb, ub = soslyapb(s, 1, factory=factory, tol=tol)
        @test log(lb) ≈ log(2.814640557) rtol=tol
        @test log(ub) ≈ log(3.980502849) rtol=tol
        lb, ub = soslyapb(s, 2, factory=factory, tol=tol)
        @test log(lb) ≈ log(3.299750624) rtol=tol
        @test log(ub) ≈ log(3.924086919) rtol=tol
        @test_throws ArgumentError sosbuildsequence(s, 1, v_0 = 2)
        for d in (1, 2)
            for p_0 in (:Primal, :Random)
                seq = sosbuildsequence(s, d, p_0=p_0)
                psw = findsmp(seq)
                @test psw !== nothing
                @test psw == smp
            end
        end
    end
    if isempty(sdp_factories)
        @test getub(s) ≈ pρub[end]
    else
        @test log(getub(s)) ≈ log(3.924086919) rtol=5e-4
        @test getsmp(s) == smp
        @test getlb(s) ≈ 3.917384715148 rtol=1e-12
    end
end

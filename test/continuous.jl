# Example 1 of
# Guglielmi, Protasov,
# "Computing Lyapunov exponents of switching systems."
# The Lyapunov Exponent was unknown
# found by SOS-24
# It is between
# (-0.7163555099121912 given by [(2,1.0786514),(1,2.467595)])
# -0.7157393916973624 given by [(2,1.0576859),(1,2.492785)] found by SOS-24
# -0.7148034620151038 given by [(2,1.02),(1,2.539)] improved from SOS-24 by hand
# and
# -0.7120361328125 given by SOS-24

@testset "Example 1" begin
    s = ContinuousSwitchedSystem([logm([1 1; -1 1]/3),logm([1 1; -1 0]/3)])
    #smp = DiscretePeriodicSwitching(s, [1, 2])
    const expected = [-0.6483688354492188, # 1
                      -0.6867463399790471, # 2
                      -0.703767519177436,  # 3
                      -0.7071168233766547, # 4
                      -0.7073838521852485, # 5
                      -0.7096574117458425, # 6
                      -0.710710268191155,  # 7
                      -0.7107560445583425, # 8
                      -0.7113358785427175, # 9
                      -0.7117890677652648, # 10
                      -0.7118256033710013, # 11
                      -0.7120087088397513] # 12
    for factory in sdp_factories
        for d in 1:12
            tol = (d < 10 ? (ismosek(factory) ? 1e-5 : 1e-4) : (ismosek(factory) ? 1e-4 : 1e-3))
            lb, ub = soslyapb(s, d, factory=factory, tol=tol)
            @test lb == -Inf
            if ismosek(factory)
                @test isapprox(ub, expected[d])
            else
                @test isapprox(ub, expected[d], rtol=2*tol)
            end
#               for k in 1:100
#                   psw = sosbuildsequence(s, d, p_0=:Random, niter=20)
#                   @test psw !== nothing
#                   @show psw.period
#                   @show psw.growthrate
#               end
        end
    end
    #@test getsmp(s) == smp
    #@test isapprox(s.lb, 3.917384715148, rtol=1e-12)
    if isempty(sdp_factories)
        @test isapprox(s.ub, Inf)
    else
        @test isapprox(s.ub, expected[12], rtol=1e-3)
        #@test_approx_eq(s.ub, -0.7107086181640625)
    end
end

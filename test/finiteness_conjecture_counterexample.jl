# See examples/Finiteness_conjecture_counterexample.ipynb

using HybridSystems
using StaticArrays
using SwitchOnSafety

@testset "Finiteness conjecture counterexample with $lp_optimizer_constructor" for lp_optimizer_constructor in lp_optimizer_constructors
    A0 = @SMatrix [1 1; 0 1]
    A1 = @SMatrix [1 0; 1 1]
    A(α) = discreteswitchedsystem([A0, α * A1])
    s1 = A(1.0)

    s = A(0.8)
    s10 = SwitchOnSafety.periodicswitching(s, [2, 1])
    psw, done, polys = invariant_polytopes(s, lp_optimizer_constructor,
        s10, tol=0.0, verbose=0)
    @test psw == s10
    @test done
    psw, done, polys = invariant_polytopes(s, lp_optimizer_constructor,
        s10, tol=1e-15, verbose=0)
    @test psw == s10
    @test !done
    psw, done, polys = invariant_polytopes(s, lp_optimizer_constructor,
        s10, tol=1e-16, verbose=0)
    @test psw == s10
    @test done
    psw, done, polys = invariant_polytopes(s, lp_optimizer_constructor,
        s10, verbose=0)
    @test psw == s10
    @test done
end

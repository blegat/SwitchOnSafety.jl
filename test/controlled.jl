# Inspired from `SetProg/test/Test/square.jl`

using LinearAlgebra, Test

using SwitchOnSafety
using Polyhedra
const Sets = SetProg.Sets
using MultivariatePolynomials
using SemialgebraicSets

function ci_square_test(
    factory, variable::SetProg.AbstractVariable, set_test, interval_test; kws...)
    □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))
    interval = polyhedron(HalfSpace([1], 1.0) ∩ HalfSpace([-1], 1))
    Δt = 0.5
    control_system = ConstrainedLinearControlDiscreteSystem(reshape([1.0], 1, 1), reshape([Δt], 1, 1), interval, interval)
    @testset "$(typeof(system))" for system in [
        ConstrainedLinearAlgebraicDiscreteSystem([1.0 Δt], [1.0 0.0], □),
        ConstrainedLinearControlDiscreteSystem([1.0 Δt; 0.0 0.0], reshape([0.0 1.0], 2, 1), □, FullSpace()),
        algebraiclift(control_system)
    ]
        set = invariant_set(system, factory, variable; verbose=0, kws...)
        set_test(set)
    end
    @testset "$(typeof(control_system))" begin
        set = invariant_set(control_system, factory, variable; verbose=0, kws...)
        interval_test(set)
    end
end

const atol = 1e-4
const rtol = 1e-4

@testset "Ellipsoid homogeneous with $factory" for factory in sdp_factories
    ci_square_test(
        factory, Ellipsoid(symmetric=true, dimension=2),
        ◯ -> begin
            @test ◯ isa Sets.Polar{Float64, Sets.EllipsoidAtOrigin{Float64}}
            @test Sets.polar(◯).Q ≈ Symmetric([1.0 -1/4; -1/4 1.0]) atol=atol rtol=rtol
        end,
        interval -> begin
            @test interval isa Sets.Polar{Float64,SetProg.Sets.EllipsoidAtOrigin{Float64}}
            @test Sets.polar(interval).Q ≈ Symmetric(ones(1, 1)) atol=atol rtol=rtol
        end,
        volume_heuristic = nth_root)
end

@testset "Ellipsoid non-homogeneous with $factory" for factory in sdp_factories
    ci_square_test(
        factory, Ellipsoid(point=SetProg.InteriorPoint([0.0, 0.0])),
        ◯ -> begin
            @test ◯ isa Sets.PerspectiveDual{Float64, Sets.Householder{Float64, Sets.ShiftedEllipsoid{Float64}, Float64}}
            z = Sets.perspective_variable(◯)
            x, y = Sets.space_variables(◯)
            ◯_dual = Sets.perspective_dual(◯)
            @test ◯_dual.p ≈ -z^2 + x^2 - x*y/2 + y^2 atol=atol rtol=rtol
            @test ◯_dual.set.Q ≈ Symmetric([1.0 -1/4; -1/4 1.0]) atol=atol rtol=rtol
            @test ◯_dual.set.b ≈ [0.0, 0.0] atol=atol rtol=rtol
            @test ◯_dual.set.β ≈ -1.0 atol=atol rtol=rtol
            @test Sets._householder(◯_dual.h) ≈ [-1.0 0.0 0.0
                                                  0.0 1.0 0.0
                                                  0.0 0.0 1.0] atol=atol rtol=rtol
        end,
        interval -> begin
            @test interval isa Sets.Translation{Sets.Polar{Float64,Sets.EllipsoidAtOrigin{Float64}},Float64,Vector{Float64}}
            @test Sets.polar(interval.set).Q ≈ Symmetric(ones(1, 1)) atol=atol rtol=rtol
            @test interval.c ≈ [0.0] atol=atol rtol=rtol
        end,
        volume_heuristic = nth_root, λ = 1.0)
end

@testset "Quadratic non-homogeneous with $factory" for factory in sdp_factories
    ci_square_test(
        factory, PolySet(degree=2, convex=true, point=SetProg.InteriorPoint([0.0, 0.0])),
        ◯ -> begin
            @test ◯ isa Sets.PerspectiveDual{Float64, Sets.Householder{Float64, Sets.ConvexPolynomialSet{Float64}, Float64}}
            z = Sets.perspective_variable(◯)
            x, y = Sets.space_variables(◯)
            ◯_dual = Sets.perspective_dual(◯)
            # The coefficient of `x*y` does not influence the volume
            # and with the values of the other parameters, it should
            # simply be in the interval [-2, -0.5].
            α = coefficient(◯_dual.p, x*y)
            @test α ≥ -2 - 2atol - rtol
            @test α ≤ -0.5 + 0.5atol + rtol
            @test ◯_dual.p ≈ -z^2 + x^2 + α*x*y + y^2 atol=atol rtol=rtol
        end,
        interval -> begin
            @test interval isa Sets.PerspectiveDual{Float64, Sets.Householder{Float64, Sets.UnknownSet{Float64}, Float64}}
            z = Sets.perspective_variable(interval)
            x, = Sets.space_variables(interval)
            interval_dual = Sets.perspective_dual(interval)
            @test interval_dual.p ≈ -z^2 + x^2 atol=atol rtol=rtol
        end,
        volume_heuristic = set -> L1_heuristic(set, [1.0, 1.0]), λ = 1.0)
end

@testset "Quartic homogeneous with $factory" for factory in sdp_factories
    ci_square_test(
        factory, PolySet(symmetric=true, degree=4, convex=true),
        ◯ -> begin
            @test ◯ isa Sets.Polar{Float64, Sets.ConvexPolynomialSublevelSetAtOrigin{Float64}}
            @test Sets.polar(◯).degree == 4
            x, y = variables(Sets.polar(◯).p)
            α = coefficient(Sets.polar(◯).p, x^3*y) / 2
            q = x^4 + 2α*x^3*y + 6x^2*y^2 + 2α*x*y^3 + y^4
            @test all(eigvals(Matrix(Sets.polar(◯).p.Q)) .≥ -atol)
            @test polynomial(Sets.polar(◯).p) ≈ q atol=atol rtol=rtol
            convexity_proof = Sets.convexity_proof(◯)
            @test convexity_proof.n == 4

            hess = 6 * [2.0, α, 2.0, α, 2.0,
                        2.0, 2.0, α, α, 2.0]
            Hess = SetProg.SumOfSquares.MultivariateMoments.SymMatrix(hess, 4)
            @test all(eigvals(Matrix(Hess)) .≥ -atol)
            @test convexity_proof.Q ≈ hess atol=atol rtol=rtol
        end,
        interval -> begin
            @test interval isa Sets.Polar{Float64, Sets.ConvexPolynomialSublevelSetAtOrigin{Float64}}
            @test Matrix(Sets.polar(interval).p.Q) ≈ ones(1, 1) atol=atol rtol=rtol
        end,
        volume_heuristic = set -> L1_heuristic(set, [1.0, 1.0]))
end

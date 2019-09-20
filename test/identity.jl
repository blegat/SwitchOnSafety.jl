# Inspired from `SetProg/test/Test/square.jl`

using LinearAlgebra, Test

using SwitchOnSafety
using Polyhedra
const Sets = SetProg.Sets

function square_test(
    factory, variable::SetProg.AbstractVariable,
    set_test; kws...)
    □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))
    for system in [
        ConstrainedContinuousIdentitySystem(2, □),
        ConstrainedDiscreteIdentitySystem(2, □)]

        set = invariant_set(system, factory, variable; verbose=0, kws...)
        set_test(set)
    end
end

const atol = 1e-5
const rtol = 1e-5

@testset "Homogeneous ellipsoid with $factory" for factory in sdp_factories
    square_test(
        factory, Ellipsoid(symmetric=true, dimension=2),
        ◯ -> begin
            @test ◯ isa Sets.Polar{Float64, Sets.EllipsoidAtOrigin{Float64}}
            @test Sets.polar(◯).Q ≈ Symmetric([1.0 0.0; 0.0 1.0]) atol=atol rtol=rtol
        end, volume_heuristic = nth_root)
end

@testset "Non-homogeneous ellipsoid with $factory" for factory in sdp_factories
    square_test(
        factory, Ellipsoid(point=SetProg.InteriorPoint([0.0, 0.0])),
        ◯ -> begin
            @test ◯ isa Sets.PerspectiveDual{Float64, Sets.Householder{Float64, Sets.ShiftedEllipsoid{Float64}, Float64}}
            z = Sets.perspective_variable(◯)
            x, y = Sets.space_variables(◯)
            ◯_dual = Sets.perspective_dual(◯)
            @test ◯_dual.p ≈ -z^2 + x^2 + y^2 atol=atol rtol=rtol
            @test Sets._householder(◯_dual.h) ≈ [-1.0 0.0 0.0
                                                  0.0 1.0 0.0
                                                  0.0 0.0 1.0] atol=atol rtol=rtol
            @test ◯_dual.set.Q ≈ Symmetric([1.0 0.0; 0.0 1.0]) atol=atol rtol=rtol
            @test ◯_dual.set.b ≈ [0.0, 0.0] atol=atol rtol=rtol
            @test ◯_dual.set.β ≈ -1.0 atol=atol rtol=rtol
        end,
        volume_heuristic = nth_root)
end

@testset "Non-homogeneous quadratic with $factory" for factory in sdp_factories
    square_test(
        factory,
        PolySet(degree=2, convex=true, point=SetProg.InteriorPoint([0.0, 0.0])),
        ◯ -> begin
            @test ◯ isa Sets.PerspectiveDual{Float64, Sets.Householder{Float64, Sets.ConvexPolynomialSet{Float64}, Float64}}
            z = Sets.perspective_variable(◯)
            x, y = Sets.space_variables(◯)
            ◯_dual = Sets.perspective_dual(◯)
            @test ◯_dual.p ≈ -z^2 + x^2 + y^2 atol=atol rtol=rtol
        end,
        volume_heuristic = set -> L1_heuristic(set, [1.0, 1.0]))
end

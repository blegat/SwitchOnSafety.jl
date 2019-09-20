# Inspired from `SetProg/test/Test/square.jl`

using LinearAlgebra, Test

using SwitchOnSafety
using Polyhedra
const Sets = SetProg.Sets
using MultivariatePolynomials

function invariant_square_test(
    factory, variable::SetProg.AbstractVariable, set_test; kws...)
    □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))
    A = [0.0 -1.0
         1.0  0.0]
    system = ConstrainedLinearDiscreteSystem(A, □)

    set = invariant_set(system, factory, variable; verbose=0, kws...)
    set_test(set)
end

const atol = 1e-4
const rtol = 1e-4

@testset "Ellipsoid homogeneous with $factory" for factory in sdp_factories
    invariant_square_test(
        factory, Ellipsoid(symmetric=true, dimension=2),
        ◯ -> begin
            @test ◯ isa Sets.Polar{Float64, Sets.EllipsoidAtOrigin{Float64}}
            @test Sets.polar(◯).Q ≈ Symmetric([1.0 0.0; 0.0 1.0]) atol=atol rtol=rtol
        end,
        volume_heuristic = nth_root)
end

@testset "Convex quadratic homogeneous with $factory" for factory in sdp_factories
    invariant_square_test(
        factory, PolySet(degree=2, convex=true, symmetric=true),
        ◯ -> begin
            @test ◯ isa Sets.Polar{Float64, Sets.ConvexPolynomialSublevelSetAtOrigin{Float64}}
            x, y = Sets.space_variables(◯)
            ◯_polar = Sets.polar(◯)
            @test ◯_polar.p ≈ x^2 + y^2 atol=atol rtol=rtol
        end,
        volume_heuristic = nth_root)
end

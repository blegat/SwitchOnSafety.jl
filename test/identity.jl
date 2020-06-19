# Inspired from `SetProg/test/Test/square.jl`

using LinearAlgebra, Test

using SwitchOnSafety
using Polyhedra
const Sets = SetProg.Sets
using MultivariatePolynomials

function square_test(
    optimizer_constructor, variable::SetProg.AbstractVariable,
    set_test; kws...)
    □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))
    for system in [
        ConstrainedContinuousIdentitySystem(2, □),
        ConstrainedDiscreteIdentitySystem(2, □)]

        set = invariant_set(system, optimizer_constructor, variable; verbose=0, kws...)
        set_test(set)
    end
end

const atol = 1e-4
const rtol = 1e-4

@testset "Homogeneous ellipsoid with $optimizer_constructor" for optimizer_constructor in sdp_factories
    square_test(
        optimizer_constructor, Ellipsoid(symmetric=true, dimension=2),
        ◯ -> begin
            @test ◯ isa Sets.Polar{Float64, Sets.Ellipsoid{Float64}}
            @test Sets.polar(◯).Q ≈ Symmetric([1.0 0.0; 0.0 1.0]) atol=atol rtol=rtol
        end, volume_heuristic = nth_root)
end

@testset "Non-homogeneous ellipsoid with $optimizer_constructor" for optimizer_constructor in sdp_factories
    square_test(
        optimizer_constructor, Ellipsoid(point=SetProg.InteriorPoint([0.0, 0.0])),
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

@testset "Non-homogeneous quadratic with $optimizer_constructor" for optimizer_constructor in sdp_factories
    square_test(
        optimizer_constructor,
        PolySet(degree=2, convex=true, point=SetProg.InteriorPoint([0.0, 0.0])),
        ◯ -> begin
            @test ◯ isa Sets.PerspectiveDual{Float64, Sets.Householder{Float64, Sets.ConvexPolynomialSet{Float64,SetProg.Sets.MonoBasis,Float64}, Float64}}
            z = Sets.perspective_variable(◯)
            x, y = Sets.space_variables(◯)
            ◯_dual = Sets.perspective_dual(◯)
            @test ◯_dual.p ≈ -z^2 + x^2 + y^2 atol=atol rtol=rtol
        end,
        volume_heuristic = set -> L1_heuristic(set, [1.0, 1.0]))
end

const quartic_inner_poly = [3.1518541833100864, -0.1617384194869734]
const quartic_inner_obj = 6.447419478140056
const quartic_inner_α = 5.6567546886722795
const quartic_inner_convexity = [12.0, 0.0, quartic_inner_α, 0.0, quartic_inner_poly[1]+2quartic_inner_poly[2],
                                 quartic_inner_α, 8.48516455194103, 0.0, 0.0, 12.0]

@testset "Quartic inner homogeneous with $optimizer_constructor" for optimizer_constructor in sdp_factories
    square_test(
        optimizer_constructor,
        PolySet(symmetric=true, degree=4, dimension=2, convex=true),
        ◯ -> begin
            @test ◯ isa Sets.Polar{Float64, Sets.ConvexPolySet{Float64,SetProg.Sets.MonoBasis,Float64}}
            ◯_polar = Sets.polar(◯)
            @test ◯_polar.degree == 4
            x, y = variables(◯_polar.p)
            @test polynomial(◯_polar.p) ≈ x^4 + quartic_inner_convexity[5]*x^2*y^2 + y^4 atol=atol rtol=rtol
            convexity_proof = Sets.convexity_proof(◯)
            @test convexity_proof.n == 4
            @test convexity_proof.Q ≈ quartic_inner_convexity atol=atol rtol=rtol
        end,
        volume_heuristic = nth_root)
end

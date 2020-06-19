# Inspired from `SetProg/test/Test/square.jl`

using LinearAlgebra, Test

using SwitchOnSafety
using Polyhedra
const Sets = SetProg.Sets
using MultivariatePolynomials

function invariant_square_test(
    optimizer_constructor, variable::SetProg.AbstractVariable, set_test; kws...)
    □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))
    A = [0.0 -1.0
         1.0  0.0]
    system = ConstrainedLinearDiscreteSystem(A, □)

    set = invariant_set(system, optimizer_constructor, variable; verbose=0, kws...)
    set_test(set)
end

const atol = 1e-4
const rtol = 1e-4

@testset "Ellipsoid homogeneous with $optimizer_constructor" for optimizer_constructor in sdp_factories
    invariant_square_test(
        optimizer_constructor, Ellipsoid(symmetric=true, dimension=2),
        ◯ -> begin
            @test ◯ isa Sets.Polar{Float64, Sets.Ellipsoid{Float64}}
            @test Sets.polar(◯).Q ≈ Symmetric([1.0 0.0; 0.0 1.0]) atol=atol rtol=rtol
        end,
        volume_heuristic = nth_root)
end

@testset "Convex quadratic homogeneous with $optimizer_constructor" for optimizer_constructor in sdp_factories
    # cholesky miss condition :: not positive definite :: info = 27 :: line 785 in sdpa_linear.cpp
    issdpa(optimizer_constructor) && continue
    invariant_square_test(
        optimizer_constructor, PolySet(degree=2, convex=true, symmetric=true),
        ◯ -> begin
            @test ◯ isa Sets.Polar{Float64, Sets.ConvexPolySet{Float64,SetProg.Sets.MonoBasis,Float64}}
            x, y = Sets.space_variables(◯)
            ◯_polar = Sets.polar(◯)
            @test ◯_polar.p ≈ x^2 + y^2 atol=atol rtol=rtol
        end,
        volume_heuristic = nth_root)
end

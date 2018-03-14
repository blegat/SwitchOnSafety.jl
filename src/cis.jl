using SemialgebraicSets
using Polyhedra

export getis

function algebraiclift(s::LinearControlDiscreteSystem)
    n = statedim(s)
    z = find(i -> iszero(sum(abs.(s.B[i,:]))), 1:n)
    @show z
    # TODO ty - 1//2y^3 + 3//1xy + 2//1yhe affine space may not be parallel to classical axis
    LinearAlgebraicDiscreteSystem(s.A[z, :], (eye(n))[z, :])
end
algebraiclift(s::DiscreteIdentitySystem) = s
algebraiclift(S::AbstractVector) = algebraiclift.(S)
algebraiclift(S::ConstantVector) = ConstantVector(algebraiclift(first(S)), length(S))
function algebraiclift(h::HybridSystem)
    HybridSystem(h.automaton, algebraiclift(h.modes), h.invariants, h.guards, algebraiclift(h.resetmaps), h.switchings)
end

function r(A::Matrix{T}, c::Vector{T}=zeros(T, size(A, 1))) where T
    [one(T) zeros(T, 1, size(A, 2))
     c A]
end

"""
    householder(x)

Householder reflection
```math
I - 2 v v^T / (v^T v)
```
It is symmetric and orthogonal.
"""
function householder(x)
    y = copy(x)
    t = LinAlg.reflector!(y)
    v = [1; y[2:end]]
    eye(length(x)) - t * v * v'
end

using DynamicPolynomials
using MultivariatePolynomials
using PolyJuMP
using SumOfSquares
using JuMP
using LightGraphs
function ATrp(p, x, A)
    B = r(A)'
    y = x[1:size(B, 2)]
    p(x => r(A)' * y)
end
function lhs(p, x, Re::ConstantVector)
    # If it is not constant, I would be comparing dummy variables of different meaning
    ATrp(p, x, first(Re).A)
end

function getp(m::Model, c, x, z::AbstractVariable)
    y = [z; x]
    n = length(x)
    #β = 1.#@variable m lowerbound=0.
    β = @variable m
    b = @variable m [1:n]
    #@constraint m b .== 0
    Q = @variable m [1:n, 1:n] Symmetric
    @constraint m y' * [β+1 b'; b Q] * y in DSOSCone()
    H = householder([1; c]) # We add 1, for z
    P = [β b'
         b Q]
    HPH = H * P * H
    p = y' * HPH * y
    ConeLyap(p, Q, b, c, H)
    #@constraint m sum(Q) == 1 # dehomogenize
    #@variable m L[1:n, 1:n]
    #@variable m λinv[1:(n-1)] >= 0
    #@SDconstraint m [Q  L
    #                 L' diagm([λinv; -1])] ⪰ 0
    #ConeLyap(x' * Q * x, Q, L, λinv)
end
function getis(s::HybridSystem{<:AbstractAutomaton, DiscreteIdentitySystem, <:LinearAlgebraicDiscreteSystem}, solver, c=map(cv->cv[1], chebyshevcenter.(s.invariants)))
    @show c
    g = s.automaton.G
    n = nv(g)
    m = SOSModel(solver=solver)
    @polyvar x[1:2] z
    y = [z; x]
    l = [getp(m, c[u], x, z) for u in 1:n]
    #X = monomials(x, 2)
    #@variable m p[1:n] Poly(X)
    #@variable m vol
    #@objective m Max vol

    @objective m Max sum(p -> trace(p.Q), l)

    λouts = Vector{Vector{JuMP.Variable}}(n)
    #λouts = Vector{Vector{Float64}}(n)
    for u in 1:n
        # Constraint 1
        N = length(out_neighbors(g, u))
        λout = @variable m [1:N] lowerbound=0
        λouts[u] = λout
        #@constraint m sum(λin) == 1
        Σ = symbol.(s.automaton, Edge.(u, out_neighbors(g, u)))
        expr = -lhs(l[u].p, y, s.resetmaps[Σ])
        for (j, v) in enumerate(out_neighbors(g, u))
            E = s.resetmaps[Σ[j]].E
            #expr -= λout[j] * p[v](x => r(E)' * x)
            newp = ATrp(l[v].p, y, E)
            expr += λout[j] * newp
        end
        @constraint m expr in DSOSCone()
        # Constraint 2
        #@SDconstraint m differentiate(p[u], x, 2) >= 0
        # Constraint 3
        for hs in ineqs(s.invariants[u])
            @constraint m l[u].p(y => [-hs.β; hs.a]) <= 0
        end
    end

    status = solve(m)

    @show getobjectivevalue(m)

    for u in 1:n
        @show getvalue.(λouts[u])
    end

    @show status
    @assert status == :Optimal
    ellipsoid.(l)
#    if uc
#        [Ellipsoid(inv(Q[u].value), c[u] + 2*reshape(C[u].value, d)) for u in vertices(g)]
#    else
#        [Ellipsoid(inv(Q[u].value), c[u]) for u in vertices(g)]
#    end
end

function getis(s::HybridSystem{<:AbstractAutomaton, DiscreteIdentitySystem, <:LinearControlDiscreteSystem}, args...)
    getis(algebraiclift(s), args...)
end

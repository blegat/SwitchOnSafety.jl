using SemialgebraicSets
using Polyhedra

export getis, fillis!, algebraiclift

function algebraiclift(s::LinearControlDiscreteSystem)
    n = statedim(s)
    z = find(i -> iszero(sum(abs.(s.B[i,:]))), 1:n)
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
# If Re is a vector not constant, I would be comparing dummy variables of different meaning
function lhs(p, x, Re::DiscreteLinearAlgebraicSystem)
    ATrp(p, x, Re.A)
end

function getp(m::Model, c, y, cone)
    n = length(y)-1
    #β = 1.#@variable m lowerbound=0.
    β = @variable m
    b = @variable m [1:n]
    #@constraint m b .== 0
    Q = @variable m [1:n, 1:n] Symmetric
    @constraint m y' * [β+1 b'; b Q] * y in cone
    H = householder([1; c]) # We add 1, for z
    P = [β b'
         b Q]
    HPH = H * P * H
    p = y' * _HPH(Q, b, β, H) * y
    vol = @variable m
    #@constraint m vol <= trace Q
    @constraint m [vol; [Q[i, j] for j in 1:n for i in 1:j]] in MOI.RootDetConeTriangle(n)
    ConeLyap(p, Q, b, β, c, H, vol)
    #@constraint m sum(Q) == 1 # dehomogenize
    #@variable m L[1:n, 1:n]
    #@variable m λinv[1:(n-1)] >= 0
    #@SDconstraint m [Q  L
    #                 L' diagm([λinv; -1])] ⪰ 0
    #ConeLyap(x' * Q * x, Q, L, λinv)
end
const DTAHAS = HybridSystem{<:AbstractAutomaton, DiscreteIdentitySystem, <:LinearAlgebraicDiscreteSystem}
function _vars(s::DTAHAS)
    @polyvar x[1:statedim(s, 1)] z
    [z; x]
end
function _p(q, N, y, l, is, ps)
    if q in N
        l[q].p
    else
        if isnull(ps[q])
            le = LiftedEllipsoid(is[q])
            ps[q] = Nullable(y' * inv(le.P) * y)
        end
        get(ps[q])
    end
end
function fillis!(is, N, s::DTAHAS, solver, c=map(cv->cv[1], chebyshevcenter.(s.invariants)); y=_vars(s), ps=fill(Nullable{polynomialtype(y, Float64)}(), length(is)), cone=DSOSCone(), λ = nothing)
    @show c
    n = nstates(s)
    m = SOSModel(solver=solver)
    l = [getp(m, c[u], y, cone) for u in N]

    @objective m Max sum(p -> p.vol, l)

    if λ === nothing
        λouts = Vector{Vector{JuMP.Variable}}(length(l))
    else
        λouts = λ
    end

    for (i, u) in enumerate(N)
        # Constraint 1
        NN = length(out_transitions(s, u))
        if λ === nothing
            λouts[i] = @variable m [1:NN] lowerbound=0
        end
        λout = λouts[i]
        for (j, t) in enumerate(out_transitions(s, u))
            σ = symbol(s.automaton, t)
            startp = lhs(l[i].p, y, s.resetmaps[σ])
            v = target(s, t)
            E = s.resetmaps[σ].E
            newp = ATrp(_p(v, N, y, l, is, ps), y, E)
            expr = λout[j] * newp - startp
            @constraint m expr in cone
        end
        # Constraint 2
        #@SDconstraint m differentiate(p[u], x, 2) >= 0
        # Constraint 3
        for hs in ineqs(s.invariants[u])
            @constraint m l[i].p(y => [-hs.β; hs.a]) <= 0
        end
    end

    JuMP.solve(m)

    @show JuMP.terminationstatus(m)
    @show JuMP.primalstatus(m)
    @show JuMP.dualstatus(m)

    @show JuMP.objectivevalue(m)

    if λ === nothing
        for i in 1:length(N)
            @show JuMP.resultvalue.(λouts[i])
        end
    end

    for (i, q) in enumerate(N)
        lv = JuMP.resultvalue(l[i])
        ps[q] = Nullable(lv.p)
        is[q] = ellipsoid(lv)
    end
end

function getis(s::DTAHAS, args...; kws...)
    nmodes = nstates(s)
    is = Vector{Ellipsoid{Float64}}(nmodes)
    fillis!(is, 1:nmodes, s, args...; kws...)
    is
end

const UnboundedControl = DiscreteLinearControlSystem{<:Any,<:Any,FullSpace}

function fillis!(is, N, s::HybridSystem{<:AbstractAutomaton, DiscreteIdentitySystem, <:HRep, FullSpace, <:UnboundedControl}, args...; kws...)
    fillis!(is, N, algebraiclift(s), args...; kws...)
end
function getis(s::HybridSystem{<:AbstractAutomaton, DiscreteIdentitySystem, <:LinearControlDiscreteSystem}, args...; kws...)
    getis(algebraiclift(s), args...; kws...)
end

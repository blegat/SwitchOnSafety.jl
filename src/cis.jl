using SemialgebraicSets
using Polyhedra
using FillArrays

export getis, fillis!, algebraiclift

function algebraiclift(s::LinearControlDiscreteSystem)
    n = statedim(s)
    z = findall(i -> iszero(sum(abs.(s.B[i,:]))), 1:n)
    # TODO ty - 1//2y^3 + 3//1xy + 2//1yhe affine space may not be parallel to classical axis
    LinearAlgebraicDiscreteSystem(s.A[z, :], Matrix(1.0I, n, n)[z, :])
end
algebraiclift(s::ConstrainedDiscreteIdentitySystem) = s
algebraiclift(S::AbstractVector) = algebraiclift.(S)
algebraiclift(S::Fill) = Fill(algebraiclift(first(S)), length(S))
function algebraiclift(h::HybridSystem)
    HybridSystem(h.automaton, algebraiclift(h.modes), algebraiclift(h.resetmaps), h.switchings)
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
    t = LinearAlgebra.reflector!(y)
    v = [1; y[2:end]]
    I - t * v * v'
end

using DynamicPolynomials
using MultivariatePolynomials
using SumOfSquares
using MathOptInterface
const MOI = MathOptInterface
using JuMP
using LightGraphs
function ATrp(p, x, A)
    B = r(A)'
    y = x[1:size(B, 2)]
    p(x => B * y)
end
# If Re is a vector not constant, I would be comparing dummy variables of different meaning
function lhs(p, x, Re::LinearAlgebraicDiscreteSystem)
    ATrp(p, x, Re.A)
end

const DTAHAS = HybridSystem{<:AbstractAutomaton, <:ConstrainedDiscreteIdentitySystem, <:LinearAlgebraicDiscreteSystem}
const DTAHCS = HybridSystem{<:AbstractAutomaton, <:ConstrainedDiscreteIdentitySystem, <:LinearControlDiscreteSystem}
function _vars(s::DTAHAS)
    @polyvar x[1:statedim(s, 1)] z
    [z; x]
end
function _p(state, N, y, l, is, ps, h)
    if state in N
        l[state].p
    else
        if ps[state] === nothing
            le = SetProg.Sets.LiftedEllipsoid(is[state])
            Pd = inv(le.P)
            H = SetProg.Sets._householder(h[state])
            HPdH = H * Pd * H
            # HPdH is not like a solution what would be obtained by solving the program
            # since the λ computed for unlifting it is maybe not one.
            # Therefore, the S-procedure's λ for the constraints will be different.
            B, b, β, λ = Bbβλ(HPdH)
            ps[state] = y' * _HPH(B/λ, b/λ, β/λ, H) * y
        end
        ps[state]
    end
end

function fillis!(is, N, s::DTAHAS, factory::JuMP.OptimizerFactory,
                 h=map(cv->InteriorPoint(cv[1]), chebyshevcenter.(stateset.(s.modes)));
                 y=_vars(s),
                 ps=fill!(Vector{Union{Nothing, polynomialtype(y, Float64)}}(undef, length(is)), nothing),
                 cone=SOSCone(),
                 λ=Dict{transitiontype(s), Float64}(),
                 enabled = 1:nstates(s),
                 detcone = MOI.RootDetConeTriangle,
                 verbose=1)
    n = nstates(s)
    model = SOSModel(factory)
    sets = @variable(model, [q in N], Ellipsoid(point=h[q]))

    @objective model Max sum(set -> nth_root(volume(set)), sets)

    λouts = Dict{transitiontype(s), JuMP.AffExpr}()

    for q in N
        @constraint(model, sets[q] ⊆ stateset(s, q))
        # Invariance constraint
        for t in out_transitions(s, q)
            λin = get(λ, t, nothing)
            target_state = target(s, t)
            if target_state in enabled
                source_set = sets[q]
                target_set = _p(target_state, N, l, is, ps, h)
                r = s.resetmaps[symbol(s, t)]
                λouts[t] = @constraint(model, r.A * source_set ⊆ r.E * target_set)
            end
        end
    end

    JuMP.optimize!(model)

    if verbose >= 1
        @show MOI.get(model, MOI.SolveTime())
        @show JuMP.termination_status(model)
        @show JuMP.primal_status(model)
        @show JuMP.dual_status(model)
        @show JuMP.objective_value(model)
    end

    if verbose >= 2
        for (t, λout) in λouts
            println("λ for $t is $(JuMP.value.(λout))")
        end
    end

    for q in N
        lv = JuMP.value(l[q])
        ps[q] = lv.p
        is[q] = ellipsoid(lv)
    end
end

function getis(s::DTAHAS, args...; kws...)
    nmodes = nstates(s)
    is = Vector{SetProg.Sets.Ellipsoid{Float64}}(undef, nmodes)
    fillis!(is, 1:nmodes, s, args...; kws...)
    is
end

function fillis!(is, N, s::DTAHCS, args...; kws...)
    fillis!(is, N, algebraiclift(s), args...; kws...)
end
function getis(s::DTAHCS, args...; kws...)
    getis(algebraiclift(s), args...; kws...)
end

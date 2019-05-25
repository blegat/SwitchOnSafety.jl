export periodicswitching, findsmp, hasrepetition
abstract type AbstractPeriodicSwitching end

# Vector{Int} for unconstrained, Vector{Edge} for constrained
duration(seq::AbstractVector) = length(seq)
#duration(seq::AbstractVector{Tuple{Int, Float64}}) = sum(map(p->p[2], seq))

adaptgrowthrate(g, len::Int) = g^(1/len)
#adaptgrowthrate(g, Δt::Float64) = log(g)/Δt
# Vector{Int} for unconstrained, Vector{Edge} for constrained
adaptgrowthrate(g, period::AbstractVector) = adaptgrowthrate(g, duration(period))
#adaptgrowthrate(g, period::AbstractVector{Tuple{Int,Float64}}) = adaptgrowthrate(g, duration(period))

mutable struct DiscretePeriodicSwitching{S<:AbstractHybridSystem, TT} <: AbstractPeriodicSwitching
    s::S
    period::Vector{TT}
    growthrate::Float64
    hashcomputed::Bool
    hash::UInt
end

function DiscretePeriodicSwitching(s, period, growthrate)
    DiscretePeriodicSwitching(s, period, growthrate, false, zero(UInt))
end

#struct ContinuousPeriodicSwitching <: AbstractPeriodicSwitching
#    s::AbstractContinuousSwitchedSystem
#    period::Vector{Tuple{Int, Float64}}
#    growthrate::Float64
#end

function Base.show(io::IO, s::AbstractPeriodicSwitching)
    #print(io, "Periodic switching of growth rate $(s.growthrate) and modes: $(symbol.(s.s, s.period))")# for the transitions $(s.period)")
    print(io, "PSW($(s.growthrate), $(symbol.(s.s, s.period)))")# for the transitions $(s.period)")
end

function periodicswitching(s::AbstractDiscreteSwitchedSystem, period::Vector, growthrate, args...)
    DiscretePeriodicSwitching(s, period, growthrate)
end

_scale(A, ::Nothing) = A
_scale(A, scaling) = A * scaling
_unscale(A, ::Nothing) = A
_unscale(A, scaling) = A / scaling

function periodicswitching(s::AbstractDiscreteSwitchedSystem, period::Vector, A::AbstractMatrix; scaling = nothing)
    lambda = ρ(A)
    growthrate = _scale(adaptgrowthrate(abs(lambda), period), scaling)
    periodicswitching(s, period, growthrate)
end

function periodicswitching(s::AbstractSwitchedSystem, period::Vector; scaling = nothing)
    A = prod(t -> _unscale(integratorfor(s, t), scaling), reverse(period))
    return periodicswitching(s, period, A, scaling = scaling)
end

# Shortcut
function periodicswitching(s::DiscreteSwitchedLinearSystem, period::Vector{Int})
    periodicswitching(s, HybridSystems.OneStateTransition.(period))
end

isperiodic(sw::HybridSystems.DiscreteSwitchingSequence) = source(sw.s, sw) == target(sw.s, sw)
function periodicswitching(s::HybridSystems.DiscreteSwitchingSequence)
    k = repetition(s.seq)
    if iszero(k)
        periodicswitching(s.s, s.seq, s.A)
    else
        periodicswitching(s.s, s.seq[1:k])
    end
end

#function periodicswitching(s::AbstractContinuousSwitchedSystem, seq, growthrate, dt)
#    # seq is a copy since it has been obtained with seq.seq[i:j]
#    mode, Δt = seq[end]
#    seq[end] = mode, dt
#    ContinuousPeriodicSwitching(s, seq, growthrate)
#end

#hasrepetition(s::AbstractPeriodicSwitching) == !iszero(repetition(s.period))

function Base.:(==)(s1::AbstractPeriodicSwitching, s2::AbstractPeriodicSwitching)
    # Shortcut
    hash(s1) == hash(s2) # assumes norepetition
#    if !(s1.s === s2.s)
#        @assert hash(s1) != hash(s2)
#        false
#    elseif !isapprox(s1.growthrate, s2.growthrate)
#        @assert hash(s1) != hash(s2)
#        false
#    else
#        p1 = s1.period
#        p2 = s2.period
#        k1 = length(p1)
#        k2 = length(p2)
#        if k1 != k2
#            k = lcm(k1, k2)
#            if k != k1
#                p1 = repmat(p1, div(k, k1))
#            end
#            if k != k2
#                p2 = repmat(p2, div(k, k2))
#            end
#        else
#            k = k1
#        end
#        for i in 1:k
#            if p1[1:i] == p2[end-i+1:end] && p1[i+1:end] == p2[1:end-i]
#                if hash(s1) != hash(s2)
#                    @show s1
#                    @show s2
#                end
#                return true
#            end
#        end
#        @assert hash(s1) != hash(s2)
#        false
#    end
end

function Base.hash(sw::DiscretePeriodicSwitching, h::UInt)
    if !sw.hashcomputed
        p = sw.period
        k = length(p)
        hs = map(i -> hash([p[end-i+1:end]; p[1:end-i]]), 1:k)
        sw.hash = hash(sort(hs))
        sw.hashcomputed = true
    end
    hash(sw.hash, h)
end

function isbetter(g1, k1, s2::AbstractPeriodicSwitching)
    g2 = s2.growthrate
    k2 = length(s2.period)
    g1 >= g2 * (1 + eps(g2)) || (g1 >= g2 * (1 - eps(g2)) && k1 < k2)
end
function isbetter(s1::AbstractPeriodicSwitching, s2::AbstractPeriodicSwitching)
    g1 = s1.growthrate
    k1 = length(s1.period)
    isbetter(g1, k1, s2)
end

periodicswitchingtype(s::AbstractDiscreteSwitchedSystem) = DiscretePeriodicSwitching{typeof(s)}
#periodicswitchingtype(s::AbstractContinuousSwitchedSystem) = ContinuousPeriodicSwitching

function bestperiod(s::AbstractDiscreteSwitchedSystem, seq::Vector, I, ::AbstractMatrix, Q::AbstractMatrix)
    periodicswitching(s, seq[I], Q), 1
end

#using Optim
#function bestperiod(s::AbstractContinuousSwitchedSystem, seq::Vector{Tuple{Int,Float64}}, I, P::AbstractMatrix, ::AbstractMatrix)
#    mode, Δt = seq[last(I)]
#    T = duration(@view seq[I]) - Δt
#    function f(dt)
#        -adaptgrowthrate(abs(ρ(integratorfor(s, (mode,dt)) * P)), T+dt)
#    end
#    res = optimize(f, min(1e-5, Δt/2), Δt, iterations = 10)
#    -Optim.minimum(res), Optim.minimizer(res)
#end

function repetition(seq)
    k = length(seq)
    for i in 1:div(k, 2)
        if iszero(k % i)
            ok = true
            for j in 2:div(k, i)
                if (@view seq[i*(j-2) .+ (1:i)]) != (@view seq[i*(j-1) .+ (1:i)])
                    ok = false
                    break
                end
            end
            if ok
                return i
            end
        end
    end
    return 0
end

"""
    findsmp(seq::HybridSystems.DiscreteSwitchingSequence)

Extract the cycle of highest growth rate in the sequence `seq`.
"""
function findsmp(seq::HybridSystems.DiscreteSwitchingSequence)
    s = seq.s
    PS = periodicswitchingtype(s)
    smp = nothing
    for i in 1:seq.len
        startNode = state(s, seq.seq[i], false)
        P = _eyet(s, seq.seq[i])
        k = 0
        for j in i:seq.len
            mode = seq.seq[j]
            k = k + nlabels(s, mode)
            Q = integratorfor(s, mode) * P
            if state(s, mode, true) == startNode
                seqij = seq.seq[i:j]
                if iszero(repetition(symbol.(s, seqij)))
                    newsmp, dt = bestperiod(s, seq.seq, i:j, P, Q)
                    notifyperiodic!(s, newsmp)
                    if smp === nothing || isbetter(newsmp, smp)
                        smp = newsmp
                    end
                end
            end
            P = Q
        end
    end

    if smp !== nothing
        updatesmp!(s, smp)
    end
    smp
end

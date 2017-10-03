import Base.start, Base.done, Base.next, Base.append!, Base.push!

abstract type AbstractSwitchingSequence end
abstract type AbstractDiscreteSwitchingSequence <: AbstractSwitchingSequence end

mutable struct DiscreteSwitchingSequence{S <: AbstractDiscreteSwitchedSystem, MT <: AbstractMatrix, TT} <: AbstractDiscreteSwitchingSequence
    s::S
    A::MT # /!\ could be dynamic or integrator
    seq::Vector{TT}
    len::Int
end

function Base.show(io::IO, s::DiscreteSwitchingSequence)
    println(io, "Discrete switching sequence of length $(s.len):")
    println(io, s.seq)
end

isperiodic(sw::DiscreteSwitchingSequence) = source(sw.s, sw) == target(sw.s, sw)
function periodicswitching(s::DiscreteSwitchingSequence)
    periodicswitching(s.s, s.seq, s.A)
end

HybridSystems.source(s::HybridSystem, seq::DiscreteSwitchingSequence) = source(s, seq.seq[1])
HybridSystems.target(s::HybridSystem, seq::DiscreteSwitchingSequence) = target(s, seq.seq[end])
# Short circuit for unconstrained system
HybridSystems.source(s::HybridSystem{OneStateAutomaton}, ::DiscreteSwitchingSequence) = 1
HybridSystems.target(s::HybridSystem{OneStateAutomaton}, ::DiscreteSwitchingSequence) = 1

#mutable struct ContinuousSwitchingSequence <: AbstractSwitchingSequence
#    s::ContinuousSwitchedSystem
#    A::AbstractMatrix # /!\ could be dynamic or integrator
#    seq::Vector{Tuple{Int, Float64}}
#    len::Int
#end

function push!(s::AbstractSwitchingSequence, el)
    if s.len < length(s.seq)
        s.seq[s.len + 1] = el
    else
        push!(s.seq, el)
    end
    s.len += 1
end

function DiscreteSwitchingSequence(s::AbstractDiscreteSwitchedSystem, A::AbstractMatrix, seq::Vector)
    DiscreteSwitchingSequence(s, A, seq, length(seq))
end
switchingsequence(s::AbstractDiscreteSwitchedSystem, A::AbstractMatrix, seq::Vector) = DiscreteSwitchingSequence(s, A, seq)
function switchingsequence(s::AbstractDiscreteSwitchedSystem, len::Int=0, v::Int=1)
    DiscreteSwitchingSequence(s, _eyes(s, v, true), Vector{transitiontype(s)}(len), 0)
end

#function ConstrainedDiscreteSwitchingSequence(s::ConstrainedDiscreteSwitchedSystem, u::Int, len=0) # TODO delete this
#    ConstrainedDiscreteSwitchingSequence(s, speye(statedim(s, u)), Vector{Edge}(len), 0)
#end

#function SwitchingSequence(s::ContinuousSwitchedSystem, len=0, v=1)
#    ContinuousSwitchingSequence(s, speye(statedim(s, v)), Vector{Tuple{Int,Float64}}(len), 0)
#end

function prepend!(s::AbstractSwitchingSequence, other::AbstractSwitchingSequence)
    s.A = s.A * other.A
    if s.len < length(s.seq)
        if s.len + other.len <= length(s.seq)
            s.seq[other.len+(1:s.len)] = s.seq[1:s.len] # Cannot use view here
            s.seq[1:other.len] = @view other.seq[1:other.len]
        else
            len = length(s.seq)
            off = max(0, len - other.len)
            append!(s.seq, @view s.seq[(off+1):s.len])
            s.seq[other.len+(1:off)] = s.seq[1:off] # Cannot use view here
            if len < other.len
                s.seq[:] = @view other.seq[1:len]
                append!(s.seq, @view other.seq[(len+1):other.len])
            else
                s.seq[1:other.len] = @view other.seq[1:other.len]
            end
        end
    else
        prepend!(s.seq, other.seq)
    end
    s.len += other.len
end

function append!(s::AbstractSwitchingSequence, other::AbstractSwitchingSequence)
    s.A = other.A * s.A
    if s.len < length(s.seq)
        if s.len + other.len <= length(s.seq)
            s.seq[(s.len+1):(s.len+other.len)] = @view other.seq[1:other.len]
        else
            off = length(s.seq) - s.len
            s.seq[(s.len+1):end] = @view other.seq[1:off]
            append!(s.seq, @view other.seq[(off+1):other.len])
        end
    else
        append!(s.seq, other.seq)
    end
    s.len += other.len
end

function measurefor(μ, s::DiscreteSwitchingSequence)
    measurefor(μ, s.s, first(s.seq))
end

# Only makes sense for discrete
function dynamicfort(s::AbstractDiscreteSwitchedSystem, sw::AbstractDiscreteSwitchingSequence)
    sw.A
end

struct SwitchingIterator{S<:AbstractDiscreteSwitchedSystem}
    s::S
    k::Int
    v0::Int
    forward::Bool
end

# Iterates over all the `forward` switching of length `k` starting at `v0`
function switchings(s::AbstractDiscreteSwitchedSystem, k::Int, v0::Int, forward=true)
    SwitchingIterator(s, k, v0, forward)
end

# nextinnode(s, v, u) = u+1
# doneinnode(s, v, u) = u >= length(s.A)
# nextoutnode(s, v, u) = u+1
# doneoutnode(s, v, u) = u >= length(s.A)

# TODO use inference instead
transitiontype(::DiscreteSwitchedLinearSystem) = Int
transitiontype(::ConstrainedDiscreteSwitchedLinearSystem) = Edge

function _next!(it, seq, modeit, modest, As, i, A)
    seq[i], modest[i] = next(modeit[i], modest[i])
    B = dynamicfort(it.s, seq[i])
    As[i] = it.forward ? B * A : A * B
    As[i], state(it.s, seq[i], it.forward)
end

function start(it::SwitchingIterator)
    k = it.k
    ET = transitiontype(it.s)
    seq = Vector{ET}(k)
    I = it.forward ? (1:k) : (k:-1:1)
    v = it.v0
    A = _eyes(it.s, v, it.forward)
    As = Vector{typeof(dynamicfort(it.s, first(transitions(it.s))))}(k)
    # modeit[i] is a list of all the possible ith mode for the (i-1)th state
    modeit = Vector{Vector{ET}}(k)
    # modest[i] is the ith state of iterator modeit[i]
    modest = Vector{Int}(k)
    for i in I
        modeit[i] = v == -1 ? ET[] : io_transitions(it.s, v, it.forward)
        modest[i] = start(modeit[i])
        if done(modeit[i], modest[i])
            v = -1
        elseif i != last(I)
            A, v = _next!(it, seq, modeit, modest, As, i, A)
        end
    end
    (modeit, modest, As, seq)
end
function done(it::SwitchingIterator, st)
    modeit, modest, _, _ = st
    I = it.forward ? (it.k:-1:1) : (1:it.k)
    for i in I
        if !done(modeit[i], modest[i])
            return false
        end
    end
    true
end
function next(it::SwitchingIterator, st)
    modeit, modest, As, seq = st
    I = it.forward ? (it.k:-1:1) : (1:it.k)
    i = -1
    for j in I
        if !done(modeit[j], modest[j])
            i = j
            break
        end
    end
    @assert i != -1
    inc = it.forward ? 1 : -1
    prev = i - inc
    A = (prev >= 1 && prev <= it.k) ? As[prev] : _eyes(it.s, 1, it.forward) # FIXME, fix we should find the right state, not put 1
    while 1 <= i <= it.k
        A, v = _next!(it, seq, modeit, modest, As, i, A)
        i += inc
        if 1 <= i <= it.k
            modeit[i] = io_transitions(it.s, v, it.forward)
            modest[i] = start(modeit[i])
        end
    end
    (switchingsequence(it.s, A, copy(seq)), (modeit, modest, As, seq))
end

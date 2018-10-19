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
    k = repetition(s.seq)
    if iszero(k)
        periodicswitching(s.s, s.seq, s.A)
    else
        periodicswitching(s.s, s.seq[1:k])
    end
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

function Base.push!(s::AbstractSwitchingSequence, el)
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

function Base.append!(s::AbstractSwitchingSequence, other::AbstractSwitchingSequence)
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
    k::Int        # length of sequence
    v0::Int       # starting mode
    forward::Bool # Is the sequence going forward or backward
end

# Iterates over all the `forward` switching of length `k` starting at `v0`
function switchings(s::AbstractDiscreteSwitchedSystem, k::Int, v0::Int, forward=true)
    SwitchingIterator(s, k, v0, forward)
end

struct SwitchingIteratorState{MT, ET}
    # modeit[i] is a list of all the possible transitions for the (i-1)th mode
    modeit::Vector{Vector{ET}}
    # modest[i] is the ith state of iterator modeit[i]
    modest::Vector{Int}
    As::Vector{MT}
    seq::Vector{ET}
    function SwitchingIteratorState{MT, ET}(k::Int) where {MT, ET}
        new{MT, ET}(Vector{Vector{ET}}(undef, k), Vector{Int}(undef, k),
                    Vector{MT}(undef, k), Vector{ET}(undef, k))
    end
end

function Base.iterate(it::SwitchingIterator)
    ET = transitiontype(it.s)
    MT = typeof(dynamicfort(it.s, first(transitions(it.s))))
    st = SwitchingIteratorState{MT, ET}(it.k)
    return complete_switching(it, st, it.forward ? 1 : it.k)
end
function Base.iterate(it::SwitchingIterator, st::SwitchingIteratorState)
    next_switching(it, st, it.forward ? it.k : 1)
end
function cur_mode(it::SwitchingIterator, st::SwitchingIteratorState, i::Int)
    j = i + (it.forward ? -1 : 1)
    if j <= 0 || j > it.k
        return it.v0
    else
        return state(it.s, st.seq[j], it.forward)
    end
end
function prev_matrix(it::SwitchingIterator, st::SwitchingIteratorState, i::Int)
    j = i + (it.forward ? -1 : 1)
    if j <= 0 || j > it.k
        return I
    else
        return st.As[j]
    end
end
function process_item_state(it::SwitchingIterator, st::SwitchingIteratorState, i::Int, item_state::Nothing)
    inc = it.forward ? 1 : -1
    return next_switching(it, st, i - inc)
end
function process_item_state(it::SwitchingIterator, st::SwitchingIteratorState, i::Int, item_state)
    st.modest[i] = item_state[2]
    st.seq[i] = item_state[1]
    B = dynamicfort(it.s, st.seq[i])
    A = prev_matrix(it, st, i)
    st.As[i] = it.forward ? B * A : A * B
    return complete_switching(it, st, i + (it.forward ? 1 : -1))
end
function next_switching(it::SwitchingIterator, st::SwitchingIteratorState, i::Int)
    inc = it.forward ? 1 : -1
    if i <= 0 || i > it.k
        return nothing
    else
        item_state = iterate(st.modeit[i], st.modest[i])
        process_item_state(it, st, i, item_state)
    end
end
function complete_switching(it::SwitchingIterator, st::SwitchingIteratorState, i::Int)
    inc = it.forward ? 1 : -1
    if i <= 0 || i > it.k
        return switchingsequence(it.s, st.As[i - inc], copy(st.seq)), st
    else
        st.modeit[i] = io_transitions(it.s, cur_mode(it, st, i), it.forward)
        item_state = iterate(st.modeit[i])
        process_item_state(it, st, i, item_state)
    end
end

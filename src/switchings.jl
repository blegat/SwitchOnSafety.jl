import Base.start, Base.done, Base.next, Base.append!, Base.push!

abstract type AbstractSwitchingSequence end
abstract type AbstractDiscreteSwitchingSequence <: AbstractSwitchingSequence end

# Vector{Int} for unconstrained, Vector{Edge} for constrained
duration(seq::AbstractVector) = length(seq)
duration(seq::AbstractVector{Tuple{Int, Float64}}) = sum(map(p->p[2], seq))

mutable struct DiscreteSwitchingSequence <: AbstractDiscreteSwitchingSequence
    s::DiscreteSwitchedSystem
    A::AbstractMatrix # /!\ could be dynamic or integrator
    seq::Vector{Int}
    len::Int
end

mutable struct ConstrainedDiscreteSwitchingSequence <: AbstractDiscreteSwitchingSequence
    s::ConstrainedDiscreteSwitchedSystem
    A::AbstractMatrix # /!\ could be dynamic or integrator
    seq::Vector{Edge}
    len::Int
end

startnode(::DiscreteSwitchingSequence) = 1
endnode(::DiscreteSwitchingSequence) = 1

startnode(seq::ConstrainedDiscreteSwitchingSequence) = startnode(seq.seq[1])
endnode(seq::ConstrainedDiscreteSwitchingSequence) = endnode(seq.seq[end])

mutable struct ContinuousSwitchingSequence <: AbstractSwitchingSequence
    s::ContinuousSwitchedSystem
    A::AbstractMatrix # /!\ could be dynamic or integrator
    seq::Vector{Tuple{Int, Float64}}
    len::Int
end

function push!(s::AbstractSwitchingSequence, el)
    if s.len < length(s.seq)
        s.seq[s.len + 1] = el
    else
        push!(s.seq, el)
    end
    s.len += 1
end

function DiscreteSwitchingSequence(s::DiscreteSwitchedSystem, A::AbstractMatrix, seq::Vector{Int})
    DiscreteSwitchingSequence(s, A, seq, length(seq))
end
SwitchingSequence(s::DiscreteSwitchedSystem, A::AbstractMatrix, seq::Vector{Int}) = DiscreteSwitchingSequence(s, A, seq)
function ConstrainedDiscreteSwitchingSequence(s::ConstrainedDiscreteSwitchedSystem, A::AbstractMatrix, seq::Vector{Edge})
    ConstrainedDiscreteSwitchingSequence(s, A, seq, length(seq))
end
SwitchingSequence(s::ConstrainedDiscreteSwitchedSystem, A::AbstractMatrix, seq::Vector{Edge}) = ConstrainedDiscreteSwitchingSequence(s, A, seq)
function SwitchingSequence(s::DiscreteSwitchedSystem, len::Int=0, v::Int=1)
    DiscreteSwitchingSequence(s, speye(dim(s)), Vector{Int}(len), 0)
end
function SwitchingSequence(s::ConstrainedDiscreteSwitchedSystem, len::Int=0, v::Int=1)
    ConstrainedDiscreteSwitchingSequence(s, speye(dim(s, v)), Vector{Edge}(len), 0)
end
function ConstrainedDiscreteSwitchingSequence(s::ConstrainedDiscreteSwitchedSystem, u::Int, len=0)
    ConstrainedDiscreteSwitchingSequence(s, speye(dim(s, u)), Vector{Edge}(len), 0)
end
function SwitchingSequence(s::ContinuousSwitchedSystem, len=0, v=1)
    ContinuousSwitchingSequence(s, speye(dim(s)), Vector{Tuple{Int,Float64}}(len), 0)
end

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
    μ[first(s.seq)]
end
function measurefor(μ, s::ConstrainedDiscreteSwitchingSequence)
  μ[s.s.eid[first(s.seq)]]
end

# Only makes sense for discrete
function dynamicfor(s::AbstractDiscreteSwitchedSystem, sw::AbstractDiscreteSwitchingSequence)
    sw.A
end

struct SwitchingIterator
    s::AbstractDiscreteSwitchedSystem
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
edgetype(::DiscreteSwitchedSystem) = Int
edgetype(::ConstrainedDiscreteSwitchedSystem) = Edge

function start(it::SwitchingIterator)
    k = it.k
    ET = edgetype(it.s)
    seq = Vector{ET}(k)
    I = it.forward ? (1:k) : (k:-1:1)
    v = it.v0
    A = speye(dim(it.s, v))
    As = Vector{eltype(it.s.A)}(k)
    # modeit[i] is a list of all the possible ith mode
    modeit = Vector{Vector{ET}}(k)
    # modest[i] is the ith state of iterator modeit[i]
    modest = Vector{Int}(k)
    for i in I
        modeit[i] = v == -1 ? ET[] : modes(it.s, v, it.forward)
        modest[i] = start(modeit[i])
        if done(modeit[i], modest[i])
            v = -1
        elseif i != last(I)
            seq[i], modest[i] = next(modeit[i], modest[i])
            As[i] = dynamicfor(it.s, seq[i]) * A
            A = As[i]
            v = state(it.s, seq[i], it.forward)
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
    inc = it.forward ? 1 : -1
    prev = i - inc
    A = (prev >= 1 && prev <= it.k) ? As[prev] : speye(dim(it.s))
    while i > 0 && i <= it.k
        seq[i], modest[i] = next(modeit[i], modest[i])
        As[i] = dynamicfor(it.s, seq[i]) * A
        A = As[i]
        v = state(it.s, seq[i], it.forward)
        i += inc
        if i > 0 && i <= it.k
            modeit[i], modest[i] = modes(it.s, v, it.forward)
        end
    end
    (SwitchingSequence(it.s, A, copy(seq)), (modeit, modest, As, seq))
end

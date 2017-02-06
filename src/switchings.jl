import Base.start, Base.done, Base.next, Base.append!, Base.push!

abstract AbstractSwitchingSequence

duration(seq::AbstractVector{Int}) = length(seq)
duration(seq::AbstractVector{Tuple{Int,Float64}}) = sum(map(p->p[2], seq))

type DiscreteSwitchingSequence <: AbstractSwitchingSequence
    s::DiscreteSwitchedSystem
    A::AbstractMatrix # /!\ could be dynamic or integrator
    seq::Vector{Int}
    len::Int
end

type ContinuousSwitchingSequence <: AbstractSwitchingSequence
    s::ContinuousSwitchedSystem
    A::AbstractMatrix # /!\ could be dynamic or integrator
    seq::Vector{Tuple{Int,Float64}}
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
function SwitchingSequence(s::DiscreteSwitchedSystem, len=0)
    DiscreteSwitchingSequence(s, speye(dim(s)), Vector{Int}(len), 0)
end
function SwitchingSequence(s::ContinuousSwitchedSystem, len=0)
    ContinuousSwitchingSequence(s, speye(dim(s)), Vector{Tuple{Int,Float64}}(len), 0)
end

function append!(s::AbstractSwitchingSequence, other::AbstractSwitchingSequence)
    s.A = other.A * s.A
    if s.len < length(s.seq)
        if s.len + other.len <= length(s.seq)
            s.seq[(s.len+1):(s.len+other.len)] = other.seq
        else
            off = length(s.seq) - s.len
            s.seq[(s.len+1):end] = @view other.seq[1:off]
            append!(s.seq, @view other.seq[(off+1):end])
        end
    else
        append!(s.seq, other.seq)
    end
    s.len += other.len
end

function measurefor(μ, s::DiscreteSwitchingSequence)
    μ[first(s.seq)]
end

# Only makes sense for discrete
function dynamicfor(s::AbstractDiscreteSwitchedSystem, sw::DiscreteSwitchingSequence)
    sw.A
end

type SwitchingIterator
    s::DiscreteSwitchedSystem
    k::Int
    v0::Int
    forward::Bool
end

function switchings(s::DiscreteSwitchedSystem, k::Int, v0::Int, forward=true)
    SwitchingIterator(s, k, v0, forward)
end

nextinnode(s, v, u) = u+1
doneinnode(s, v, u) = u >= length(s.A)
nextoutnode(s, v, u) = u+1
doneoutnode(s, v, u) = u >= length(s.A)

function start(it::SwitchingIterator)
    k = it.k
    seq = Vector{Int}(k)
    I = it.forward ? (1:k) : (k:-1:1)
    v = it.v0
    A = speye(dim(it.s))
    As = Vector{eltype(it.s.A)}(k)
    modeit = Vector{Any}(k)
    modest = Vector{Any}(k)
    for i in I
        modeit[i] = v == -1 ? (1:0) : modes(it.s, v, it.forward)
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
    (DiscreteSwitchingSequence(it.s, A, copy(seq)), (modeit, modest, As, seq))
end

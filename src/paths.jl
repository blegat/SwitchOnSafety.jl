type SwitchingSequence
    s::SwitchedSystem
    A::AbstractMatrix
    seq::Vector{Int}
end

type SwitchingIterator
  s::SwitchedSystem
  k::Int
  v0::Int
  rev::Bool
end

function switchings(s::SwitchedSystem, k::Int)
    SwitchingIterator(s, k)
end

nextinnode(s, v, u) = u+1
doneinnode(s, v, u) = u >= length(s.A)
nextoutnode(s, v, u) = u+1
doneoutnode(s, v, u) = u >= length(s.A)

function start(it::SwitchingIterator)
    seq = SwitchingSequence(it.s, speye(dim(it.s)), zeros(Int, it.k))
    if it.rev
        off = 1
        I = (k-1):-1:1
        f = nextinnode
        seq.seq[k] = it.v0
    else
        off = -1
        I = 2:k
        f = nextoutnode
        seq.seq[1] = it.v0
    end
    for i in I
        seq.seq[i] = f(it.s, seq.seq[i+off], 0)
    end
    seq
end
function done(it::SwitchingIterator, state)
    if it.rev
        off = 1
        I = (k-1):-1:1
    else
        off = -1
        I = 2:k
    end
    for i in I
    end
end
function next(it::SwitchingIterator, state)
end

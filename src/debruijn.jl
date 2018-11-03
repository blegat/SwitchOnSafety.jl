export mdependentlift

"""
    mdependentlift(s::AbstractDiscreteSwitchedSystem, m, forward=true)

Returns the `m`-path-dependent lift of the hybrid system `s` [LD06]. The
corresponding automaton is the debruijn graph of dimension `m`. If `forward` is
`false` then the debruijn graph works backward [ELD14].

[LD06] Lee, J.-W. & Dullerud, G. E.
Uniform stabilization of discrete-time switched and Markovian jump linear systems
Automatica, Elsevier, 2006, 42, 205-21

[ELD14] Essick, R. and Lee, J.-W. and Dullerud, G. E.
Control of linear switched systems with receding horizon modal information
IEEE Transactions on Automatic Control, 2014, 59, 2340-2352

[PEDJ16] Philippe, M. and Essick, R. and Dullerud, G. E. and Jungers, R. M.
Stability of discrete-time switching systems with constrained switching sequences
Automatica, 2016, 72, 242-250
"""
function mdependentlift(s::AbstractDiscreteSwitchedSystem, n::Integer, forward=true)
    curstid = 0
    idmap = Dict{Int, Int}()
    statelabels = String[]
    for st in states(s)
        for sw in switchings(s, n, st, forward)
            id = stateid(s, sw.seq, 1:n)
            if !haskey(idmap, id)
                curstid += 1
                idmap[id] = curstid
                push!(statelabels, join(string.(sw.seq)))
            end
        end
    end
    G = LightAutomaton(curstid)
    for st in states(s)
        for sw in switchings(s, n+1, st, forward)
            seq = sw.seq
            src = stateid(s, seq, 1:n)
            dst = stateid(s, seq, 1 .+ (1:n))
            src_id = idmap[src]
            dst_id = idmap[dst]
            # The label is automatically the same so we just have to check
            # whether there is any transition between `src_id` and `dst_id`
            if !has_transition(G, src_id, dst_id)
                σ = forward ? symbol(s, seq[n+1]) : symbol(s, seq[1])
                add_transition!(G, src_id, dst_id, σ)
            end
        end
    end
    s = discreteswitchedsystem(map(rm -> rm.A, s.resetmaps), G)
    s.ext[:statelabels] = statelabels
    s
end

function stateid(s, seq, I)
    m = length(s.resetmaps)
    id = 0
    for i in I
        id = id * m + symbol(s, seq[i])-1
    end
    id + 1
end

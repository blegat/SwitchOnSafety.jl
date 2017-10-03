using GraphPlot
export aplot

function aplot(s, A::LightAutomaton)
    el = map(t -> symbol(s, t), transitions(A))
    if haskey(s.ext, :statelabels)
        gplot(A.G, nodelabel=s.ext[:statelabels], edgelabel=el)
    else
        gplot(A.G, edgelabel=el)
    end
end
function aplot(s, osa::OneStateAutomaton)
    A = LightAutomaton(1)
    for σ in 1:ntransitions(osa)
        add_transition!(A, 1, 1, σ)
    end
    aplot(s, A)
end
function aplot(s)
    aplot(s, s.automaton)
end

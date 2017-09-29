export findsmp
abstract type AbstractPeriodicSwitching end

adaptgrowthrate(g, len::Int) = g^(1/len)
adaptgrowthrate(g, Δt::Float64) = log(g)/Δt
# Vector{Int} for unconstrained, Vector{Edge} for constrained
adaptgrowthrate(g, period::AbstractVector) = adaptgrowthrate(g, duration(period))
adaptgrowthrate(g, period::AbstractVector{Tuple{Int,Float64}}) = adaptgrowthrate(g, duration(period))

function (::Type{T})(s::AbstractSwitchedSystem, period::Vector) where T<:AbstractPeriodicSwitching
    A = speye(statedim(s, state(s, first(period), false)))
    for edge in period
        A = integratorfor(s, edge) * A
    end
    lambda = ρ(A)
    growthrate = adaptgrowthrate(abs(lambda), period)
    T(s, period, growthrate)
end

function (==)(s1::AbstractPeriodicSwitching, s2::AbstractPeriodicSwitching)
    if !(s1.s === s2.s)
        false
    elseif !isapprox(s1.growthrate, s2.growthrate)
        false
    else
        p1 = s1.period
        p2 = s2.period
        k1 = length(p1)
        k2 = length(p2)
        if k1 != k2
            k = lcm(k1, k2)
            if k != k1
                p1 = repmat(p1, div(k, k1))
            end
            if k != k2
                p2 = repmat(p2, div(k, k2))
            end
        else
            k = k1
        end
        for i in 1:k
            if p1[1:i] == p2[end-i+1:end] && p1[i+1:end] == p2[1:end-i]
                return true
            end
        end
        false
    end
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

function updatesmp!(s::AbstractSwitchedSystem, smp::AbstractPeriodicSwitching)
    updatelb!(s, smp.growthrate)
    if isnull(s.smp) || isbetter(smp, get(s.smp))
        s.smp = smp
    end
    smp
end

function getsmp(s::AbstractSwitchedSystem)
    if isnull(s.smp)
        error("No smp found")
    end
    get(s.smp)
end

function findsmp(seq)
    s = seq.s
    PS = periodicswitchingtype(s)
    smp = Nullable{PS}()
    for i in 1:seq.len
        startNode = state(s, seq.seq[i], false)
        P = speye(statedim(s, startNode))
        k = 0
        for j in i:seq.len
            mode = seq.seq[j]
            k = k + nlabels(s, mode)
            Q = integratorfor(s, mode) * P
            if state(s, mode, true) == startNode
                growthrate, dt = bestperiod(s, seq.seq, i:j, P, Q)
                if isnull(smp) || isbetter(growthrate, length(i:j), get(smp))
                    smp = Nullable{PS}(periodicswitching(s, seq.seq[i:j], growthrate, dt))
                end
            end
            P = Q
        end
    end

    if !isnull(smp)
        updatesmp!(s, get(smp))
    end
    smp
end

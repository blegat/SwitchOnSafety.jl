abstract AbstractPeriodicSwitching

function (::Type{T}){T<:AbstractPeriodicSwitching}(s::AbstractSwitchedSystem, period::Vector)
    A = speye(dim(s))
    for mode in period
        A = matrixforward(s, mode) * A
    end
    lambda = ρ(A)
    growthrate = abs(lambda)^(1/length(period))
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
        throw(InvalidStateException("No smp found"))
    end
    get(s.smp)
end

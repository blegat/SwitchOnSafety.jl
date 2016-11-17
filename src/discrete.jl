export DiscreteSwitchedSystem, DiscretePeriodicSwitching, getsmp

abstract AbstractDiscreteSwitchedSystem <: AbstractSwitchedSystem

type DiscretePeriodicSwitching <: AbstractPeriodicSwitching
    s::AbstractDiscreteSwitchedSystem
    period::Vector{Int}
    growthrate::Float64
end
function DiscretePeriodicSwitching(s::AbstractDiscreteSwitchedSystem, period::Vector{Int})
    A = speye(dim(s))
    for mode in period
        A = matrixfor(s, mode) * A
    end
    lambda = ρ(A)
    growthrate = abs(lambda)^(1/length(period))
    DiscretePeriodicSwitching(s, period, growthrate)
end

function (==)(s1::DiscretePeriodicSwitching, s2::DiscretePeriodicSwitching)
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

function isbetter(g1, k1, s2::DiscretePeriodicSwitching)
    g2 = s2.growthrate
    k2 = length(s2.period)
    g1 >= g2 * (1 + eps(g2)) || (g1 >= g2 * (1 - eps(g2)) && k1 < k2)
end

function isbetter(s1::DiscretePeriodicSwitching, s2::DiscretePeriodicSwitching)
    g1 = s1.growthrate
    k1 = length(s1.period)
    isbetter(g1, k1, s2)
end

type DiscreteSwitchedSystem <: AbstractDiscreteSwitchedSystem
    A::Vector
    lb::Float64
    ub::Float64
    # There will typically only be lyapunov for small d so a dictionary would be overkill
    lyaps::Vector{Nullable{Lyapunov}}
    smp::Nullable{DiscretePeriodicSwitching}
    function DiscreteSwitchedSystem(A::Vector)
        if isempty(A)
            error("Needs at least one matrix in the system")
        end
        n = size(A[1], 1)
        for M in A
            if !isa(M, AbstractMatrix)
                error("One of the matrices is of invalid type: $(typeof(M))")
            end
            if LinAlg.checksquare(M) != n
                error("The matrices should all have the same dimensions")
            end
        end
        new(A, 0, Inf, Nullable{Lyapunov}[], nothing)
    end
end

function updatesmp!(s::AbstractDiscreteSwitchedSystem, smp::DiscretePeriodicSwitching)
    updatelb!(s, smp.growthrate)
    if isnull(s.smp) || isbetter(smp, get(s.smp))
        s.smp = smp
    end
    smp
end

function getsmp(s::AbstractDiscreteSwitchedSystem)
    if isnull(s.smp)
        throw(InvalidStateException("No smp found"))
    end
    get(s.smp)
end

function quicklb(s::DiscreteSwitchedSystem)
    qlb = maximum(map(ρ, s.A))
    updatelb!(s, qlb)
end

function quickub(s::DiscreteSwitchedSystem)
    qub = minimum(map(p -> maximum(map(A->norm(A,p), s.A)), [1, 2, Inf]))
    updateub!(s, qub)
end

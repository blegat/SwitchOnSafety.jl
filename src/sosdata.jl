export getlb, getub, sosdata, getsmp, hassmp

abstract type AbstractSOSData end
mutable struct SOSData{S} <: AbstractSOSData
    x::Vector{PolyVar{true}}
    lb::Float64
    ub::Float64
    # There will typically only be lyapunov for small d so a dictionary would be overkill
    lyaps::Vector{Nullable{Lyapunov}}
    smp::Nullable{DiscretePeriodicSwitching{S}}
end
function newsosdata(s::DiscreteSwitchedLinearSystem)
    n = statedim(s, 1)
    @polyvar x[1:n]
    SOSData{typeof(s)}(x, 0, Inf, Nullable{Lyapunov}[], nothing)
end

import LightGraphs
mutable struct ConstrainedSOSData{S} <: AbstractSOSData
    x::Vector{Vector{PolyVar{true}}}
    eid::Dict{LightGraphs.Edge, Int}
    lb::Float64
    ub::Float64
    lyaps::Vector{Nullable{Lyapunov}}
    smp::Nullable{DiscretePeriodicSwitching{S}}
end
function newsosdata(s::ConstrainedDiscreteSwitchedLinearSystem)
    eid = Dict{LightGraphs.Edge, Int}()
    neid = 0
    for t in transitions(s)
        neid += 1
        eid[t] = neid
    end
    y = Vector{Vector{PolyVar{true}}}(nstates(s))
    for st in states(s)
        @polyvar x[1:statedim(s, st)]
        y[st] = x
    end
    ConstrainedSOSData{typeof(s)}(y, eid, 0, Inf, Nullable{Vector{Lyapunov}}[], nothing)
end

const sosdatakey = :SwitchOnSafetyData

function sosdata(s::AbstractSwitchedSystem)
    if !haskey(s.ext, sosdatakey)
        s.ext[sosdatakey] = newsosdata(s)
    end
    s.ext[sosdatakey]
end

#states(s::AbstractSwitchedSystem) = 1:1
#nstates(s::AbstractSwitchedSystem) = 1
#statedim(s::AbstractSwitchedSystem, i::Int) = s.n[i]
MultivariatePolynomials.variables(s::SOSData, state::Int) = s.x
MultivariatePolynomials.variables(s::ConstrainedSOSData, state::Int) = s.x[state]
MultivariatePolynomials.variables(s::AbstractSwitchedSystem, state::Int) = variables(sosdata(s), state)

getlyaps(s::AbstractSOSData) = s.lyaps
getlyaps(s::AbstractSwitchedSystem) = getlyaps(sosdata(s))

getlb(s::AbstractSOSData) = s.lb
getlb(s::AbstractSwitchedSystem) = getlb(sosdata(s))
function updatelb!(s::AbstractSOSData, lb, smp=nothing)
    s.lb = max(s.lb, lb)
    lb
end
updatelb!(s::AbstractSwitchedSystem, lb, args...) = updatelb!(sosdata(s), lb, args...)

getub(s::AbstractSOSData) = s.ub
getub(s::AbstractSwitchedSystem) = getub(sosdata(s))
function updateub!(s::AbstractSOSData, ub)
    s.ub = min(s.ub, ub)
    ub
end
updateub!(s::AbstractSwitchedSystem, ub) = updateub!(sosdata(s), ub)

function updateb!(s::AbstractSwitchedSystem, lb, ub)
    updatelb!(s, lb), updateub!(s, ub)
end

function updatesmp!(s::AbstractSOSData, smp::AbstractPeriodicSwitching)
    updatelb!(s, smp.growthrate)
    if isnull(s.smp) || isbetter(smp, get(s.smp))
        s.smp = smp
    end
    smp
end
updatesmp!(s::AbstractSwitchedSystem, smp::AbstractPeriodicSwitching) = updatesmp!(sosdata(s), smp)

hassmp(s::AbstractSOSData) = !isnull(s.smp)
hassmp(s::AbstractSwitchedSystem) = hassmp(sosdata(s))
function getsmp(s::AbstractSOSData)
    if !hassmp(s)
        error("No smp found")
    end
    get(s.smp)
end
getsmp(s::AbstractSwitchedSystem) = getsmp(sosdata(s))

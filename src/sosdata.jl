export getlb, getub, sosdata, getsmp, hassmp, unstable_periodic_switchings

mutable struct Lyapunov{PT <: HybridSystems.StateProperty,
                        DT <: HybridSystems.TransitionProperty}
    d::Int
    soslb::Float64
    dual::DT
    sosub::Float64
    primal::PT
end

mutable struct SOSData{S, TT,
                       XT <: HybridSystems.StateProperty,
                       PT <: HybridSystems.StateProperty,
                       DT <: HybridSystems.TransitionProperty}
    x::XT
    lb::Float64
    ub::Float64
    # There will typically only be lyapunov for small d so a dictionary would be overkill
    lyaps::Vector{Union{Nothing, Lyapunov{PT, DT}}}
    smp::Union{Nothing, DiscretePeriodicSwitching{S, TT}}
    # Unstable periodic switchings
    upsw::Union{Nothing, Set{DiscretePeriodicSwitching{S, TT}}}
end
function SOSData{S, TT, XT, PT, DT}(s::S) where {S, TT,
                                                 XT <: HybridSystems.StateProperty,
                                                 PT <: HybridSystems.StateProperty,
                                                 DT <: HybridSystems.TransitionProperty}
    y = HybridSystems.state_property(s, Vector{PolyVar{true}})::XT
    for st in states(s)
        @polyvar x[1:statedim(s, st)]
        y[st] = x
    end
    lyaps = Union{Nothing, Lyapunov{PT, DT}}[]
    SOSData{S, TT, XT, PT, DT}(y, 0, Inf, lyaps, nothing, nothing)
end

const PolynomialLyapunov{T} = SetProg.Sets.PolynomialSublevelSetAtOrigin{T}
const MeasureLyapunov{T} = MultivariateMoments.MomentMatrix{T, DynamicPolynomials.Monomial{true}, DynamicPolynomials.MonomialVector{true}}

function SOSData(s::AbstractDiscreteSwitchedSystem)
    S = typeof(s)
    TT = transitiontype(s)
    XT = HybridSystems.state_property_type(S, Vector{PolyVar{true}})
    PT = HybridSystems.state_property_type(S, PolynomialLyapunov{Float64})
    DT = HybridSystems.transition_property_type(S, MeasureLyapunov{Float64})
    SOSData{S, TT, XT, PT, DT}(s)
end

const sosdatakey = :SwitchOnSafetyData

function sosdata(s::AbstractSwitchedSystem)
    if !haskey(s.ext, sosdatakey)
        s.ext[sosdatakey] = SOSData(s)
    end
    s.ext[sosdatakey]
end

#states(s::AbstractSwitchedSystem) = 1:1
#nstates(s::AbstractSwitchedSystem) = 1
#statedim(s::AbstractSwitchedSystem, i::Int) = s.n[i]
# using MP.variables triggers the invalid range update Julia v0.6 bug
MultivariatePolynomials.variables(s::SOSData, state::Int) = s.x[state]
MultivariatePolynomials.variables(s::AbstractSwitchedSystem, state::Int) = variables(sosdata(s), state)

getlyaps(s::SOSData) = s.lyaps

getlb(s::SOSData) = s.lb
function updatelb!(s::SOSData, lb, smp=nothing)
    s.lb = max(s.lb, lb)
    lb
end

getub(s::SOSData) = s.ub
function updateub!(s::SOSData, ub)
    s.ub = min(s.ub, ub)
    ub
end

function updateb!(s::AbstractSwitchedSystem, lb, ub)
    updatelb!(s, lb), updateub!(s, ub)
end

function updatesmp!(s::SOSData, smp::AbstractPeriodicSwitching)
    updatelb!(s, smp.growthrate)
    if s.smp === nothing || isbetter(smp, s.smp)
        s.smp = smp
    end
    smp
end

function notifyperiodic!(s::SOSData, psw::AbstractPeriodicSwitching, tol=sqrt(eps(Float64)))
    if s.upsw !== nothing && psw.growthrate+tol >= 1
        if !(psw in s.upsw)
            push!(s.upsw, psw)
        end
    end
end

function compute_unstable_periodic_switchings!(s::SOSData{S, TT}) where {S, TT}
    s.upsw = Set{DiscretePeriodicSwitching{S, TT}}()
end
unstable_periodic_switchings(s::SOSData) = s.upsw

hassmp(s::SOSData) = s.smp !== nothing
function getsmp(s::SOSData)
    if !hassmp(s)
        error("No smp found")
    end
    return s.smp
end

for f in (:getsmp, :hassmp, :unstable_periodic_switchings, :notifyperiodic!, :updatesmp!, :updateub!, :getub, :updatelb!, :getlyaps, :getlb)
    @eval begin
        $f(s::AbstractSwitchedSystem, args...) = $f(sosdata(s), args...)
    end
end

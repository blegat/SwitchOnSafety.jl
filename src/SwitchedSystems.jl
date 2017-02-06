module SwitchedSystems

using MultivariatePolynomials

import Base.==

export AbstractSwitchedSystem
export dim, ρ, quicklb, quickub, quickb

function ρ(A::AbstractMatrix)
    maximum(abs.(eigvals(A)))
end

abstract AbstractSwitchedSystem

type Lyapunov
    d::Int
    soslb::Float64
    dual::Vector{Measure{true, Float64}}
    sosub::Float64
    primal::Polynomial{true, Float64}
end

function dim(s::AbstractSwitchedSystem)
    size(s.A[1], 1)
end

function nlabels(s::AbstractSwitchedSystem, mode)
    1
end

function modes(s::AbstractSwitchedSystem, v, forward=true)
    1:length(s.A)
end

function dynamicfor(s::AbstractSwitchedSystem, mode::Int)
    s.A[mode]
end

function state(s::AbstractSwitchedSystem, mode, forward=true)
    1
end

function updatelb!(s::AbstractSwitchedSystem, lb, smp=nothing)
    s.lb = max(s.lb, lb)
    lb
end

function updateub!(s::AbstractSwitchedSystem, ub)
    s.ub = min(s.ub, ub)
    ub
end

function updateb!(s::AbstractSwitchedSystem, lb, ub)
    updatelb!(s, lb), updateub!(s, ub)
end

function quickb(s::AbstractSwitchedSystem)
    (quicklb(s), quickub(s))
end

include("periodic.jl")
include("discrete.jl")
include("continuous.jl")
include("switchings.jl")
include("veronese.jl")
include("pradius.jl")
include("sos.jl")

end # module

module SwitchedSystems

using MultivariatePolynomials
using LightGraphs

import Base.==

export AbstractSwitchedSystem
export dim, ρ, quicklb, quickub, quickb

# eigvals is not defined for SparseMatrixCSC
ρ(A::AbstractSparseMatrix) = ρ(full(A))
function ρ(A::AbstractMatrix)
    maximum(abs.(eigvals(A)))
end

abstract AbstractSwitchedSystem

type Lyapunov
    d::Int
    soslb::Float64
    dual::Vector{Measure{true, Float64}}
    sosub::Float64
    primal::Vector{Polynomial{true, Float64}}
end

nnodes(s::AbstractSwitchedSystem) = 1
dim(s::AbstractSwitchedSystem) = dim(s, 1)
dim(s::AbstractSwitchedSystem, i::Int) = s.n[i]
MultivariatePolynomials.vars(s::AbstractSwitchedSystem, i::Int) = s.x
startnode(i::Int) = 1
endnode(i::Int) = 1

nlabels(s::AbstractSwitchedSystem, edge::Int) = 1

modes(s::AbstractSwitchedSystem, v::Int, forward=true) = 1:length(s.A)

dynamicfor(s::AbstractSwitchedSystem, mode::Int) = s.A[mode]

state(s::AbstractSwitchedSystem, edge::Int, forward=true) = 1

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
include("constrained.jl")
include("continuous.jl")
include("switchings.jl")
include("veronese.jl")
include("kronecker.jl")
include("pradius.jl")
include("sos.jl")

end # module

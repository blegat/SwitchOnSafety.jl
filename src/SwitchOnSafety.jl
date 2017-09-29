module SwitchOnSafety

using DynamicPolynomials
using MultivariatePolynomials
using MultivariateMoments

using LightGraphs

import Base.==

export AbstractSwitchedSystem
export statedim, ρ, quicklb, quickub, quickb

# eigvals is not defined for SparseMatrixCSC
ρ(A::AbstractSparseMatrix) = ρ(full(A))
function ρ(A::AbstractMatrix)
    maximum(abs.(eigvals(A)))
end

abstract type AbstractSwitchedSystem end

mutable struct Lyapunov
    d::Int
    soslb::Float64
    dual::Vector # TODO measure type
    sosub::Float64
    primal::Vector # TODO polynomial type
end

states(s::AbstractSwitchedSystem) = 1:1
nstates(s::AbstractSwitchedSystem) = 1
statedim(s::AbstractSwitchedSystem, i::Int) = s.n[i]
MultivariatePolynomials.variables(s::AbstractSwitchedSystem, i::Int) = s.x
startnode(i::Int) = 1
endnode(i::Int) = 1

nlabels(s::AbstractSwitchedSystem, edge) = 1

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
include("continuous.jl")
include("switchings.jl")
include("veronese.jl")
include("kronecker.jl")
include("pradius.jl")
include("sos.jl")
include("sosseq.jl")
include("sosext.jl")

end # module

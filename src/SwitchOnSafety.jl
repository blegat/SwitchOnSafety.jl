module SwitchOnSafety

using DynamicPolynomials
using MultivariatePolynomials
using MultivariateMoments

using HybridSystems

import Base.==

export ρ, quicklb, quickub, quickb

function ρ(A::Matrix)
    maximum(abs.(eigvals(A)))
end
# eigvals is not defined for SparseMatrixCSC
ρ(A::AbstractSparseMatrix) = ρ(full(A))
# eigvals is not defined for SMatrix in StaticArrays for non-Hermitian
# I don't want to add StaticArrays to the REQUIRE file though so I
# define it for AbstractMatrix.
ρ(A::AbstractMatrix) = ρ(Matrix(A))

#abstract type AbstractSwitchedSystem end

mutable struct Lyapunov
    d::Int
    soslb::Float64
    dual::Vector # TODO measure type
    sosub::Float64
    primal::Vector # TODO polynomial type
end

const AbstractDiscreteSwitchedSystem = Union{DiscreteSwitchedLinearSystem, ConstrainedDiscreteSwitchedLinearSystem}
const AbstractSwitchedSystem = AbstractDiscreteSwitchedSystem
integratorfor(s::AbstractDiscreteSwitchedSystem, t) = dynamicfort(s, t)
#integratorfor(s::AbstractContinuousSwitchedSystem, mode::Tuple{Int,Float64}) = expm(dynamicfor(s, mode[1]) * mode[2])

ρA(s::ConstrainedDiscreteSwitchedLinearSystem) = ρ(adjacency_matrix(s.automaton.G))

state(s::AbstractSwitchedSystem, t, forward=true) = forward ? target(s, t) : source(s, t)

nlabels(s::AbstractSwitchedSystem, t) = 1

#modes(s::AbstractSwitchedSystem, v::Int, forward=true) = 1:length(s.A)

dynamicforσ(s::AbstractDiscreteSwitchedSystem, σ) = s.resetmaps[σ].A
dynamicfort(s::AbstractDiscreteSwitchedSystem, t) = dynamicforσ(s, symbol(s, t))

io_transitions(s, st, forward::Bool) = forward ? out_transitions(s, st) : in_transitions(s, st)

function _eyet(s, t)
    # transpose needed for rectangle system
    eye(integratorfor(s, t)')
end
function _eyes(s, st, forward)
    _eyet(s, first(io_transitions(s, st, forward)))
end

function quickb(s::AbstractSwitchedSystem)
    (quicklb(s), quickub(s))
end

ρA(s::DiscreteSwitchedLinearSystem) = ntransitions(s)

function quicklb(s::DiscreteSwitchedLinearSystem)
    qlb = maximum(ρ.(dynamicfort.(s, transitions(s))))
    updatelb!(s, qlb)
end

function quickub(s::AbstractDiscreteSwitchedSystem)
    qub = minimum(map(p -> maximum(norm.(dynamicfort.(s, transitions(s)), p)), [1, 2, Inf]))
    updateub!(s, qub)
end

include("periodic.jl")
include("sosdata.jl")
include("switchings.jl")
include("veronese.jl")
include("kronecker.jl")
include("pradius.jl")
include("sos.jl")
include("sosseq.jl")
include("sosext.jl")

include("debruijn.jl")

#include("discrete.jl")
#include("continuous.jl")

end # module

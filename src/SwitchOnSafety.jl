module SwitchOnSafety

using LinearAlgebra, SparseArrays, Statistics

using DynamicPolynomials
using MultivariatePolynomials
const MP = MultivariatePolynomials
using SemialgebraicSets
using MultivariateMoments

import Reexport
Reexport.@reexport using SetProg
Reexport.@reexport using MathematicalSystems
Reexport.@reexport using HybridSystems

export ρ, quicklb, quickub, quickb

struct AB
    A::Matrix{Float64}
    B::Matrix{Float64}
end

Base.:*(ab1::AB, ab2::AB) = AB(ab1.A * ab2.A, [ab1.A * ab2.B ab1.B])

function LinearAlgebra.opnorm(ab::AB, p::Real=2)
    if p != 2
        throw(ArgumentError("`opnorm(::AB, $p)` not supported, only `2` is supported."))
    end
    C = nullspace(ab.B')
    return opnorm(C' * ab.A, 2)
end

function ρ(A::Matrix)
    if isempty(A)
        throw(ArgumentError("Cannot compute the spectral radius of 0x0 matrix"))
    end
    maximum(abs.(eigvals(A)))
end

function ρ(A::Matrix, B::Matrix)
    n = size(A, 1)
    powerA = I
    Cspace = zeros(n, 0)
    for i in (1: n)
        Cspace = [Cspace (powerA * B)]
        powerA *= A
    end 
    C = nullspace(Cspace')
    E = (C) * (C)'
    A22 = (C)' * A * C
    if size(C, 2) == 0
        return 0
    end
    ρ(A22)
end

ρ(ab::AB) = ρ(ab.A, ab.B)

# eigvals is not defined for SparseMatrixCSC
# eigvals is not defined for SMatrix in StaticArrays for non-Hermitian
ρ(A::AbstractMatrix) = ρ(Matrix(A))

# TODO move to StaticArrays
using StaticArrays
function ρ(A::SMatrix{2, 2, <:Real})
    trace = tr(A)
    deter = det(A)
    # Need to solve λ^2 - trace * λ + deter
    Δ = trace^2 - 4deter
    if Δ >= 0
        sΔ = √Δ
        return max(abs(trace + sΔ), abs(trace - sΔ)) / 2
    else
        return √(trace^2 + abs(Δ)) / 2
    end
end

#abstract type AbstractSwitchedSystem end

include("scaled.jl")

const DiscreteSwitchedLinearControlSystem = HybridSystem{OneStateAutomaton, <:ContinuousIdentitySystem, <:LinearControlMap, AutonomousSwitching}
const AbstractDiscreteSwitchedLinearSystem = Union{DiscreteSwitchedLinearSystem, ConstrainedDiscreteSwitchedLinearSystem}
const AbstractDiscreteSwitchedSystem = Union{AbstractDiscreteSwitchedLinearSystem, DiscreteSwitchedLinearControlSystem}
const AbstractSwitchedSystem = AbstractDiscreteSwitchedSystem
integratorfor(s::AbstractDiscreteSwitchedSystem, t) = dynamicfort(s, t)
#integratorfor(s::AbstractContinuousSwitchedSystem, mode::Tuple{Int,Float64}) = expm(dynamicfor(s, mode[1]) * mode[2])

ρA(s::ConstrainedDiscreteSwitchedLinearSystem) = ρ(adjacency_matrix(s.automaton.G))

state(s::AbstractSwitchedSystem, t, forward=true) = forward ? target(s, t) : source(s, t)

nlabels(s::AbstractSwitchedSystem, t) = 1

#modes(s::AbstractSwitchedSystem, v::Int, forward=true) = 1:length(s.A)

function measurefor(μ, s::HybridSystems.DiscreteSwitchingSequence)
    μ[first(s.seq)]
end

# Only makes sense for discrete
function dynamicfort(s::AbstractDiscreteSwitchedSystem, sw::HybridSystems.DiscreteSwitchingSequence)
    sw.A
end

function identity_for_mode(s::AbstractDiscreteSwitchedLinearSystem, mode)
    return Matrix(LinearAlgebra.I, HybridSystems.statedim(s, mode), HybridSystems.statedim(s, mode))
end
function identity_for_mode(s::DiscreteSwitchedLinearControlSystem, mode)
    return AB(Matrix(1.0I, HybridSystems.statedim(s, mode), HybridSystems.statedim(s, mode)), zeros(HybridSystems.statedim(s, mode), 0))
end
dynamicforσ(s::AbstractDiscreteSwitchedLinearSystem, σ) = s.resetmaps[σ].A
dynamicforσ(s::DiscreteSwitchedLinearControlSystem, σ) = AB(s.resetmaps[σ].A, s.resetmaps[σ].B)
dynamicfort(s::AbstractDiscreteSwitchedSystem, t) = dynamicforσ(s, symbol(s, t))
dynamicfort(s::AbstractDiscreteSwitchedSystem, t, γ) = dynamicfort(s, t) / γ

io_transitions(s, st, forward::Bool) = forward ? out_transitions(s, st) : in_transitions(s, st)

function _eyet(s, t)
    # transpose needed for rectangle system
    n, m = size(integratorfor(s, t)')
    Matrix(1.0I, n, m)
end

include("periodic.jl")
include("quick.jl")
include("gripenberg.jl")

include("polytopes.jl")

include("sosdata.jl")
include("veronese.jl")
include("kronecker.jl")
include("pradius.jl")
include("sos.jl")
include("sosseq.jl")
include("sosext.jl")

include("debruijn.jl")

include("plot.jl")

include("invariant.jl")

end # module

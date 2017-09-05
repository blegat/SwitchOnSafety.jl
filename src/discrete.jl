export DiscreteSwitchedSystem, DiscretePeriodicSwitching, getsmp

abstract type AbstractDiscreteSwitchedSystem <: AbstractSwitchedSystem end

struct DiscretePeriodicSwitching <: AbstractPeriodicSwitching
    s::AbstractDiscreteSwitchedSystem # Cannot say DiscreteSwitchedSystem as it would be a circular type declaration https://github.com/JuliaLang/julia/issues/269
    period::Vector{Int}
    growthrate::Float64
end

integratorfor(s::AbstractDiscreteSwitchedSystem, edge) = dynamicfor(s, edge)

function bestperiod(s::AbstractDiscreteSwitchedSystem, seq::Vector, I, ::AbstractMatrix, Q::AbstractMatrix)
    adaptgrowthrate(abs(ρ(Q)), @view seq[I]), 1
end

mutable struct DiscreteSwitchedSystem <: AbstractDiscreteSwitchedSystem
    A::Vector
    n::Int
    x::Vector{PolyVar{true}}
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
        @polyvar x[1:n]
        for M in A
            if !isa(M, AbstractMatrix)
                error("One of the matrices is of invalid type: $(typeof(M))")
            end
            if LinAlg.checksquare(M) != n
                error("The matrices should all have the same dimensions")
            end
        end
        new(A, n, x, 0, Inf, Nullable{Lyapunov}[], nothing)
    end
end

ρA(s::DiscreteSwitchedSystem) = length(s.A)

function quicklb(s::DiscreteSwitchedSystem)
    qlb = maximum(map(ρ, s.A))
    updatelb!(s, qlb)
end

function quickub(s::AbstractDiscreteSwitchedSystem)
    qub = minimum(map(p -> maximum(map(A->norm(A, p), s.A)), [1, 2, Inf]))
    updateub!(s, qub)
end

periodicswitchingtype(s::DiscreteSwitchedSystem) = DiscretePeriodicSwitching

function periodicswitching(s::DiscreteSwitchedSystem, seq::Vector{Int})
    gr = adaptgrowthrate(ρ(prod(reverse(s.A[seq]))), seq)
    periodicswitching(s, seq, gr)
end
function periodicswitching(s::DiscreteSwitchedSystem, seq::Vector{Int}, growthrate, args...)
    DiscretePeriodicSwitching(s, seq, growthrate)
end

include("constrained.jl")

export DiscreteSwitchedSystem, DiscretePeriodicSwitching, getsmp

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
discreteswitchedsystem(A) = DiscreteSwitchedSystem(A)

periodicswitchingtype(s::DiscreteSwitchedSystem) = DiscretePeriodicSwitching

function periodicswitching(s::DiscreteSwitchedSystem, seq::Vector{Int})
    gr = adaptgrowthrate(ρ(prod(reverse(s.A[seq]))), seq)
    periodicswitching(s, seq, gr)
end
function periodicswitching(s::DiscreteSwitchedSystem, seq::Vector{Int}, growthrate, args...)
    DiscretePeriodicSwitching(s, seq, growthrate)
end

include("constrained.jl")

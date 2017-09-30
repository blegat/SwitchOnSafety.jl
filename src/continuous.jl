export ContinuousSwitchedSystem, ContinuousPeriodicSwitching, getsmp

abstract type AbstractContinuousSwitchedSystem <: AbstractSwitchedSystem end

mutable struct ContinuousSwitchedSystem <: AbstractContinuousSwitchedSystem
    A::Vector
    n::Int
    x::Vector{PolyVar{true}}
    lb::Float64
    ub::Float64
    # There will typically only be lyapunov for small d so a dictionary would be overkill
    lyaps::Vector{Nullable{Lyapunov}}
    smp::Nullable{ContinuousPeriodicSwitching}
    function ContinuousSwitchedSystem(A::Vector)
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
        new(A, n, x, -Inf, Inf, Nullable{Lyapunov}[], nothing)
    end
end

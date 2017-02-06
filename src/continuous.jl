using Optim

export ContinuousSwitchedSystem, ContinuousPeriodicSwitching, getsmp

abstract AbstractContinuousSwitchedSystem <: AbstractSwitchedSystem

type ContinuousPeriodicSwitching <: AbstractPeriodicSwitching
    s::AbstractContinuousSwitchedSystem
    period::Vector{Tuple{Int, Float64}}
    growthrate::Float64
end
integratorfor(s::AbstractContinuousSwitchedSystem, mode::Tuple{Int,Float64}) = expm(dynamicfor(s, mode[1]) * mode[2])


function bestperiod(s::AbstractContinuousSwitchedSystem, seq::Vector{Tuple{Int,Float64}}, I, P::AbstractMatrix, ::AbstractMatrix)
    mode, Δt = seq[last(I)]
    T = duration(@view seq[I]) - Δt
    function f(dt)
        -adaptgrowthrate(abs(ρ(integratorfor(s, (mode,dt)) * P)), T+dt)
    end
    res = optimize(f, min(1e-5, Δt/2), Δt, iterations = 10)
    -Optim.minimum(res), Optim.minimizer(res)
end

type ContinuousSwitchedSystem <: AbstractContinuousSwitchedSystem
    A::Vector
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
        for M in A
            if !isa(M, AbstractMatrix)
                error("One of the matrices is of invalid type: $(typeof(M))")
            end
            if LinAlg.checksquare(M) != n
                error("The matrices should all have the same dimensions")
            end
        end
        new(A, -Inf, Inf, Nullable{Lyapunov}[], nothing)
    end
end

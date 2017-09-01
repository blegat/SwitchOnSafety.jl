using Optim

export ContinuousSwitchedSystem, ContinuousPeriodicSwitching, getsmp

abstract type AbstractContinuousSwitchedSystem <: AbstractSwitchedSystem end

struct ContinuousPeriodicSwitching <: AbstractPeriodicSwitching
    s::AbstractContinuousSwitchedSystem
    period::Vector{Tuple{Int, Float64}}
    growthrate::Float64
end
integratorfor(s::AbstractContinuousSwitchedSystem, mode::Tuple{Int,Float64}) = expm(dynamicfor(s, mode[1]) * mode[2])

function nullsmp(s::AbstractContinuousSwitchedSystem)
    Nullable{ContinuousPeriodicSwitching}()
end
function buildsmp(s::AbstractContinuousSwitchedSystem, seq, growthrate, dt)
    # seq is a copy since it has been obtained with seq.seq[i:j]
    mode, Δt = seq[end]
    seq[end] = mode, dt
    Nullable{ContinuousPeriodicSwitching}(ContinuousPeriodicSwitching(s, seq, growthrate))
end
function bestperiod(s::AbstractContinuousSwitchedSystem, seq::Vector{Tuple{Int,Float64}}, I, P::AbstractMatrix, ::AbstractMatrix)
    mode, Δt = seq[last(I)]
    T = duration(@view seq[I]) - Δt
    function f(dt)
        -adaptgrowthrate(abs(ρ(integratorfor(s, (mode,dt)) * P)), T+dt)
    end
    res = optimize(f, min(1e-5, Δt/2), Δt, iterations = 10)
    -Optim.minimum(res), Optim.minimizer(res)
end

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

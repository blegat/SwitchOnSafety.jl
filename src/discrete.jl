export DiscreteSwitchedSystem, DiscretePeriodicSwitching, getsmp

abstract AbstractDiscreteSwitchedSystem <: AbstractSwitchedSystem

type DiscretePeriodicSwitching <: AbstractPeriodicSwitching
    s::AbstractDiscreteSwitchedSystem
    period::Vector{Int}
    growthrate::Float64
end
matrixforward(s::AbstractDiscreteSwitchedSystem, mode::Int) = matrixfor(s, mode)

type DiscreteSwitchedSystem <: AbstractDiscreteSwitchedSystem
    A::Vector
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
        for M in A
            if !isa(M, AbstractMatrix)
                error("One of the matrices is of invalid type: $(typeof(M))")
            end
            if LinAlg.checksquare(M) != n
                error("The matrices should all have the same dimensions")
            end
        end
        new(A, 0, Inf, Nullable{Lyapunov}[], nothing)
    end
end

function quicklb(s::DiscreteSwitchedSystem)
    qlb = maximum(map(ρ, s.A))
    updatelb!(s, qlb)
end

function quickub(s::DiscreteSwitchedSystem)
    qub = minimum(map(p -> maximum(map(A->norm(A,p), s.A)), [1, 2, Inf]))
    updateub!(s, qub)
end

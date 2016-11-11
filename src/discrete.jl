export SwitchedSystem

type SwitchedSystem <: AbstractSwitchedSystem
    A::Vector
    lb::Float64
    ub::Float64
    # There will typically only be lyapunov for small d so a dictionary would be overkill
    lyaps::Vector{Nullable{Lyapunov}}
    function SwitchedSystem(A::Vector)
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
        new(A, 0, Inf, Nullable{Lyapunov}[])
    end
end

function quicklb(s::SwitchedSystem)
    qlb = maximum(map(ρ, s.A))
    updatelb!(s, qlb)
end

function quickub(s::SwitchedSystem)
    qub = minimum(map(p -> maximum(map(A->norm(A,p), s.A)), [1, 2, Inf]))
    updateub!(s, qub)
end

module SwitchedSystems

export SwitchedSystem
export dim, ρ, quicklb, quickub, quickb

function ρ(A::AbstractMatrix)
  maximum(abs.(eigvals(A)))
end

type SwitchedSystem
  A::Vector
  lb::Float64
  ub::Float64
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
    new(A, 0, Inf)
  end
end

function dim(s::SwitchedSystem)
  size(s.A[1], 1)
end

function updatelb!(s, lb)
  s.lb = max(s.lb, lb)
end

function updateub!(s, ub)
  s.ub = min(s.ub, ub)
end

function quicklb(s::SwitchedSystem)
  qlb = maximum(map(ρ, s.A))
  updatelb!(s, qlb)
  qlb
end

function quickub(s::SwitchedSystem)
  qub = minimum(map(p -> maximum(map(A->norm(A,p), s.A)), [1, 2, Inf]))
  updateub!(s, qub)
  qub
end

function quickb(s::SwitchedSystem)
  (quicklb(s), quickub(s))
end

include("veronese.jl")
include("pradius.jl")
include("sos.jl")

end # module

module SwitchedSystems

export SwitchedSystem
export ρ, quicklb, quickub, quickb

type SwitchedSystem
  A::Vector
  lb::Float64
  ub::Float64
  function SwitchedSystem(A::Vector)
    new(A, 0, Inf)
  end
end

function ρ(A::AbstractMatrix)
  maximum(abs.(eigvals(A)))
end

function quicklb(s::SwitchedSystem)
  qlb = maximum(map(ρ, s.A))
  s.lb = max(s.lb, qlb)
  qlb
end

function quickub(s::SwitchedSystem)
  qub = minimum(map(p -> maximum(map(A->norm(A,p), s.A)), [1, 2, Inf]))
  s.ub = min(s.ub, qub)
  qub
end

function quickb(s::SwitchedSystem)
  (quicklb(s), quickub(s))
end

include("veronese.jl")

end # module

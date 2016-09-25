module SwitchedSystems

export SwitchedSystem
export ρ, quicklb, quickub, quickb

type SwitchedSystem
  A::Vector
  jsrlb::Float64
  jsrub::Float64
  function SwitchedSystem(A::Vector)
    new(A, 0, Inf)
  end
end

function ρ(A::AbstractMatrix)
  maximum(abs.(eigvals(A)))
end

function quicklb(s::SwitchedSystem)
  qlb = maximum(map(ρ, s.A))
  s.jsrlb = max(s.jsrlb, qlb)
  qlb
end

function quickub(s::SwitchedSystem)
  qub = minimum(map(p -> maximum(map(A->norm(A,p), s.A)), [1, 2, Inf]))
  s.jsrub = min(s.jsrub, qub)
  qub
end

function quickb(s::SwitchedSystem)
  (quicklb(s), quickub(s))
end

end # module

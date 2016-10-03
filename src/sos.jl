using JuMP
using SumOfSquares
using MathProgBase

export sosb, soscert

function soscert(s::SwitchedSystem, d, γ; solver=MathProgBase.defaultSDPsolver)
  n = dim(s)
  @polyvar x[1:n]
  m = JuMP.Model(solver=solver)
  Z = monomials(x, 2*d)
  @SOSvariable m p Z
  @SOSconstraint m p >= sum(x.^(2*d))
  for A in s.A
    @SOSconstraint m p(A*x, x) <= γ^(2*d) * p
  end
  status = solve(m)
  if status == :Optimal
    status, getvalue(p), nothing
  elseif status == :Infeasible
    status, nothing, nothing
  else
    error("Solver returned with status : $status")
  end
end

function sosb(s::SwitchedSystem, d::Integer; solver=MathProgBase.defaultSDPsolver, tol=1e-5)
  # The SOS ub is greater than the JSR hence also greater than any of its lower bound
  soslb = s.lb
  (lb, sosub) = pradiusb(s, 2*d)
  updatelb!(s, lb)
  while sosub - soslb > tol
    mid = (soslb + sosub) / 2
    status, primal, dual = soscert(s, d, mid; solver=solver)
    if status == :Optimal
      sosub = mid
    elseif status == :Infeasible
      soslb = mid
    end
  end
  ub = sosub
  n = dim(s)
  lb = sosub / min(length(s.A), binomial(n+d-1, d))^(1/(2*d))
  updatelb!(s, lb)
  updateub!(s, ub)
  lb, ub
end

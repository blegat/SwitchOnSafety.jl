function sosduallyap(s::SwitchedSystem, d, γ; solver=MathProgBase.defaultSDPsolver)
  n = dim(s)
  @polyvar x[1:n]
  model = JuMP.Model(solver=solver)
  Z = monomials(x, 2*d)
  @SOSvariable model p Z
  @SOSconstraint model p >= sum(x.^(2*d))
  cons = [@SOSconstraint model p(A*x, x) >= γ^(2*d) * p for A in s.A]
  status = solve(model)
  if status == :Optimal
    status, getvalue(p), nothing
  elseif status == :Infeasible
    status, nothing, nothing #[getdual(c) for c in cons]
  else
    error("Solver returned with status : $status")
  end
end


function sosduallyapb(s::SwitchedSystem, d::Integer; solver=MathProgBase.defaultSDPsolver, tol=1e-5)
  # The SOS ub is greater than the JSR hence also greater than any of its lower bound
  soslb = 0
  sosub = 10
  while sosub - soslb > tol
      @show sosub
      @show soslb
    mid = (soslb + sosub) / 2
    status, primal, dual = sosduallyap(s, d, mid; solver=solver)
    if status == :Optimal
      soslb = mid
    elseif status == :Infeasible
      sosub = mid
    end
  end
  sosub
# n = dim(s)
# lb = sosub / min(length(s.A), binomial(n+d-1, d))^(1/(2*d))
# updatelb!(s, lb)
# updateub!(s, ub)
# lb, ub
end

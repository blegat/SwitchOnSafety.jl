function sosduallyap(s::SwitchedSystem, d::Integer, γ; solver=MathProgBase.defaultSDPsolver)
  n = dim(s)
  m = length(s.A)
  @polyvar x[1:n]
  model = JuMP.Model(solver=solver)
  Z = monomials(x, d)
  p = Vector{MatPolynomial{JuMP.Variable}}(m)
  for i in 1:m
      @SOSvariable model tmp >= 0 Z
      p[i] = tmp
  end
  lhs = 0
  rhs = 0
  for i in 1:m
      lhs += p[i](inv(s.A[i])*x, x)
      rhs += p[i]
  end
  @SOSconstraint model lhs >= γ^(2*d) * rhs
  lhs = 0
  for i in 1:n
      y = zeros(n)
      y[i] = 1
      lhs += sum([q(y, x) for q in p])
  end
  @show lhs
  @SOSconstraint model lhs >= 1
  status = solve(model)
  if status == :Optimal
    status, map(getvalue, p), nothing
  elseif status == :Infeasible
    status, nothing, nothing
  else
    error("Solver returned with status : $status")
  end
end

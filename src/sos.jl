using JuMP
using SumOfSquares
using MathProgBase

export soslyap, soslyapb, sosduallyap, sosduallyapb

function setlyap!(s, lyap::Lyapunov)
  d = lyap.d
  if length(s.lyaps) < d
    sizehint!(s.lyaps, d)
    while length(s.lyaps) < d
      push!(s.lyaps, nothing)
    end
  end
  s.lyaps[d] = lyap
end

function getlyap(s::SwitchedSystem, d::Int; solver=MathProgBase.defaultSDPsolver, tol=1e-5)
  if d > length(s.lyaps) || isnull(s.lyaps[d])
    soslyapb(s, d; solver=solver, tol=1e-5, cached=true)
  end
  get(s.lyaps[d])
end


function soslyap(s::SwitchedSystem, d, γ; solver=MathProgBase.defaultSDPsolver)
  n = dim(s)
  @polyvar x[1:n]
  model = JuMP.Model(solver=solver)
  Z = monomials(x, 2*d)
  @SOSvariable model p Z
  @SOSconstraint model p >= sum(x.^(2*d))
  cons = [@SOSconstraint model p(A*x, x) <= γ^(2*d) * p for A in s.A]
  status = solve(model)
  if status == :Optimal
    status, getvalue(p), nothing
  elseif status == :Infeasible
    status, nothing, [getdual(c) for c in cons]
  else
    error("Solver returned with status : $status")
  end
end

function soslyapb(s::SwitchedSystem, d::Integer; solver=MathProgBase.defaultSDPsolver, tol=1e-5, cached=true)
  # The SOS ub is greater than the JSR hence also greater than any of its lower bound
  soslb = s.lb
  (lb, sosub) = pradiusb(s, 2*d)
  while sosub - soslb > tol
    mid = (soslb + sosub) / 2
    status, primal, dual = soslyap(s, d, mid; solver=solver)
    if status == :Optimal
      sosub = mid
    elseif status == :Infeasible
      soslb = mid
    end
  end
  if cached
    setlyap!(s, Lyapunov(d, soslb, dual, sosub, primal))
  end
  ub = sosub
  n = dim(s)
  lb = sosub / min(length(s.A), binomial(n+d-1, d))^(1/(2*d))
  updateb!(s, lb, ub)
end

function

function sosbuildsequence(s::SwitchedSystem, d::Integer; v_0=:Random, p_0=:Random, solver=MathProgBase.defaultSDPsolver, tol=1e-5)
  lyap = getlyap(s, d, solver=solver, tol=tol)
  if p_0 == :Primal
    p_0 = lyap.primal
  elseif p_0 == :Random
    @polyvar x[1:n]
    Z = monomials(x, d)
    p_0 = randsos(Z, :Gram)
  end # otherwise p_0 is assumed to be an sos polynomial given by the user
  p_0 = VecPolynomial(p_0)
  if v_0 == :Random
    curnode = randi(1)
  else
    if !isa(v_0, Integer) || v_0 < 1 || v_0 > 1
      throw(ArgumentError("Invalid v_0=$v_0"))
    end
    curnode = v_0
  end

  p_k = p_0
  n = dim(s)
  prod = eye(n)
  seq = zeros(1,niter)

  for iter = 1:l:(niter-l+1)
      best = 0
      best_e = -1
      npaths = 0

      i = l
      pathnodes = zeros(Int,l+1)
      pathnodes[l+1] = curnode
      pathedges = zeros(Int,l)
      pathprods = cell(1,l+1)
      pathprods[l+1] = eye(n)
      while i <= l
          if i == 0
              npaths = npaths + 1
              pathprod = pathprods[1]
              cur = dot(lyap.dual[pathedges[1]], p_k(pathprod*x, x))
              if cur > best
                  best = cur
                  best_e = pathedges
              end
              i = i + 1
          else
              pathedges[i] = pathedges[i] + 1
              if pathedges[i] > nEdges
                  pathedges[i] = 0
                  i = i + 1
              elseif edges(pathedges[i], 2) == pathnodes[i+1]
                  pathnodes[i] = edges(pathedges[i],1)
                  pathprods[i] = pathprods[i+1] * getmat(s, pathedges[i])
                  i = i - 1
              end
          end
      end
      if npaths == 0
          error("$curnode does not have any incoming path of length $l.")
      end
      if best_e == -1
          error("Oops, this should not happend, please report this bug.")
      end
      curnode = edges(best_e(1), 1)
      seq((niter-iter-l+2):(niter-iter+1)) = best_e
      Medges = Sys.getMatProdEdge(best_e)
      p_k = (Medges') * p_k * Medges
      prod = prod * Medges
  end

  timeBuilding = toc(initBuilding) ;

  initExtracting = tic ;

  lbd = 0 ;
  for i = 1:numel(seq)
      P = eye(Sys.getMatSize()) ;
      startNode = edges(seq(i),1) ;
      k = 0 ;
      for j = i:numel(seq)
          e = seq(j) ;
          k = k + Sys.getNLabels(e) ;
          P = Sys.getMatProdEdge(e) * P ;
          if edges(e,2) == startNode
              lambda = eigs(P,1) ;
              newLbd = abs(lambda)^(1/k) ;
              if newLbd > lbd * (1 + sqrt(eps))
                  smp = seq(i:j) ;
                  if opts.verbose > 0
                      fprintf('Better S.M.P. detected with lb = %f > %f\n', newLbd^(1/Sys.degree), lbd^(1/Sys.degree)) ;
                      fprintf('Path (edges) : %s\n', num2str(smp)) ;
                  end
                  lbd = newLbd ;
              end
          end
      end
  end

  lb = lbd^(1/Sys.degree) ;

  timeExtracting = toc(initExtracting) ;

  info = struct() ;
  info.timeBuilding = timeBuilding ;
  info.timeExtracting = timeExtracting ;

end

#function sosbuildsequence(s::SwitchedSystem, d::Integer; solver=MathProgBase.defaultSDPsolver, tol=1e-5)
#end

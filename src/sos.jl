using JuMP
using PolyJuMP
using SumOfSquares
using MathProgBase

export soslyap, soslyapb, sosbuildsequence

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
    @polyvariable model p Z
    @polyconstraint model p >= sum(x.^(2*d))
    cons = [@polyconstraint model p(A*x, x) <= γ^(2*d) * p for A in s.A]
    # I suppress the warning "Not solved to optimality, status: Infeasible"
    status = solve(model, suppress_warnings=true)
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
    primal = dual = nothing
    while sosub - soslb > tol
        mid = (soslb + sosub) / 2
        status, curprimal, curdual = soslyap(s, d, mid; solver=solver)
        if status == :Optimal
            primal = curprimal
            sosub = mid
        elseif status == :Infeasible
            dual = curdual
            soslb = mid
        end
    end
    if cached
        if primal === nothing
            status, primal, _ = soslyap(s, d, solub; solver=solver)
            @assert status == :Optimal
            @assert !(primal === nothing)
        end
        if dual === nothing
            status, _, dual = soslyap(s, d, soslb; solver=solver)
            @assert status == :Infeasible
            @assert !(dual === nothing)
        end
        setlyap!(s, Lyapunov(d, soslb, dual, sosub, primal))
    end
    ub = sosub
    n = dim(s)
    lb = sosub / min(length(s.A), binomial(n+d-1, d))^(1/(2*d))
    updateb!(s, lb, ub)
end

function sosbuildsequence(s::SwitchedSystem, d::Integer; v_0=:Random, p_0=:Random, l=1, niter=42, solver=MathProgBase.defaultSDPsolver, tol=1e-5)
    lyap = getlyap(s, d, solver=solver, tol=tol)
    if p_0 == :Primal
        p_0 = lyap.primal
        x = vars(p_0)
    elseif p_0 == :Random
        @polyvar x[1:n]
        Z = monomials(x, d)
        p_0 = randsos(Z, :Gram)
    end # otherwise p_0 is assumed to be an sos polynomial given by the user
    p_0 = VecPolynomial(p_0)
    if v_0 == :Random
        curstate = rand(1:1)
    else
        if !isa(v_0, Integer) || v_0 < 1 || v_0 > 1
            throw(ArgumentError("Invalid v_0=$v_0"))
        end
        curstate = v_0
    end

    p_k = p_0
    n = dim(s)
    prod = speye(n)
    seq = zeros(Int, niter)

    for iter = 1:l:(niter-l+1)
        best = 0
        best_seq = nothing
        nswitchings = 0
        for cur_seq in switchings(s, l, curstate, false)
            nswitchings = nswitchings + 1
            cur = dot(lyap.dual[first(cur_seq.seq)], p_k(cur_seq.A * x, x))
            @assert cur_seq.A == s.A[cur_seq.seq[1]]
            if cur > best
                best = cur
                best_seq = cur_seq
            end
        end

        if nswitchings == 0
            error("$curstate does not have any incoming path of length $l.")
        end
        if best_seq == nothing
            error("Oops, this should not happend, please report this bug.")
        end
        curstate = state(s, best_seq.seq[1], false)
        seq[(niter-iter-l+2):(niter-iter+1)] = best_seq.seq
        p_k = p_k(best_seq.A * x, x)
        prod = prod * best_seq.A
    end

    smp = Nullable{PeriodicSwitching}()
    for i = 1:length(seq)
        P = speye(n)
        startNode = state(s, seq[i], false)
        k = 0
        for j = i:length(seq)
            mode = seq[j]
            k = k + nlabels(s, mode)
            P = matrixfor(s, mode) * P
            if state(s, mode, true) == startNode
                lambda = ρ(P)
                growthrate = abs(lambda)^(1/k)
                if isnull(smp) || isbetter(growthrate, length(i:j), get(smp))
                    smp = Nullable{PeriodicSwitching}(PeriodicSwitching(s, seq[i:j], growthrate))
                end
            end
        end
    end

    if !isnull(smp)
        updatesmp!(s, get(smp))
    end
    smp
end

#function sosbuildsequence(s::SwitchedSystem, d::Integer; solver=MathProgBase.defaultSDPsolver, tol=1e-5)
#end

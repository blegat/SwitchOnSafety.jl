using JuMP
using PolyJuMP
using SumOfSquares
using MathProgBase

export soslyap, soslyapb, sosbuildsequence

# Storing the Lyapunov
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
function getlyap(s::AbstractSwitchedSystem, d::Int; solver=MathProgBase.defaultSDPsolver, tol=1e-5)
    if d > length(s.lyaps) || isnull(s.lyaps[d])
        soslyapb(s, d; solver=solver, tol=1e-5, cached=true)
    end
    get(s.lyaps[d])
end

# Building the Lyapunov constraints
function soslyapforward(s::AbstractDiscreteSwitchedSystem, p, A)
    x = vars(p)
    p(A*x, x)
end
function soslyapforward(s::AbstractContinuousSwitchedSystem, p, A)
    x = vars(p)
    dot(differentiate(p, x), A*x)
end
soslyapscaling(s::AbstractDiscreteSwitchedSystem, γ, d) = γ^(2*d)
soslyapscaling(s::AbstractContinuousSwitchedSystem, γ, d) = 2*d*γ
function soslyapconstraint(s::AbstractSwitchedSystem, model::JuMP.Model, p, A, d, γ)
    @polyconstraint model soslyapforward(s, p, A) <= soslyapscaling(s, γ, d) * p
end
function soslyapconstraints(s::AbstractSwitchedSystem, model::JuMP.Model, p, d, γ)
    [soslyapconstraint(s, model, p, A, d, γ) for A in s.A]
end

# Solving the Lyapunov problem
function soslyap(s::AbstractSwitchedSystem, d, γ; solver=MathProgBase.defaultSDPsolver)
    n = dim(s)
    @polyvar x[1:n]
    model = JuMP.Model(solver=solver)
    Z = monomials(x, 2*d)
    @polyvariable model p Z
    @polyconstraint model p >= sum(x.^(2*d))
    cons = soslyapconstraints(s, model, p, d, γ)
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

function getsoslyapinitub(s::AbstractDiscreteSwitchedSystem, d::Integer)
    _, sosub = pradiusb(s, 2*d)
    sosub
end
function getsoslyapinitub(s::AbstractContinuousSwitchedSystem, d::Integer)
    Inf
end

function increaselb(s::AbstractDiscreteSwitchedSystem, lb, step)
    lb *= step
end

soschecktol(soslb, sosub, tol) = sosub - soslb > tol
soschecktol(s::AbstractDiscreteSwitchedSystem, soslb, sosub, tol) = soschecktol(log(soslb), log(sosub), tol)
soschecktol(s::AbstractContinuousSwitchedSystem, soslb, sosub, tol) = soschecktol(soslb, sosub, tol)

function sosmid(soslb, sosub, step)
    if isfinite(soslb) && isfinite(sosub)
        mid = (soslb + sosub) / 2
    elseif isfinite(soslb)
        mid = soslb + step
    elseif isfinite(sosub)
        mid = sosub - step
    else
        mid = 0
    end
end
sosmid(s::AbstractDiscreteSwitchedSystem, soslb, sosub, step) = exp(sosmid(log(soslb), log(sosub), step))
sosmid(s::AbstractContinuousSwitchedSystem, soslb, sosub, step) = sosmid(soslb, sosub, step)

function soslb2lb(s::AbstractDiscreteSwitchedSystem, sosub, d)
    n = dim(s)
    η = min(length(s.A), binomial(n+d-1, d))
    sosub / η^(1/(2*d))
end
soslb2lb(s::AbstractContinuousSwitchedSystem, soslb, d) = -Inf

# Obtaining bounds with Lyapunov
function soslyapb(s::AbstractSwitchedSystem, d::Integer; solver=MathProgBase.defaultSDPsolver, tol=1e-5, step=1, cached=true)
    # The SOS ub is greater than the JSR hence also greater than any of its lower bound
    soslb = s.lb
    sosub = getsoslyapinitub(s, d)
    primal = dual = nothing
    while soschecktol(s, soslb, sosub, tol)
        mid = sosmid(s, soslb, sosub, step)
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
    lb = soslb2lb(s, soslb, d)
    updateb!(s, lb, ub)
end

function candidates(s::AbstractDiscreteSwitchedSystem, l, curstate)
    switchings(s, l, curstate, false)
end
function candidates(s::AbstractContinuousSwitchedSystem, l, curstate)
    modes(s, curstate, false)
end
measurefor(μs, dyn::Int) = μs[dyn]

function best_dynamic(s::AbstractSwitchedSystem, μs, p::VecPolynomial, l, curstate)
    best = 0
    best_dyn = nothing
    ncandidates = 0
    for dyn in candidates(s, l, curstate)
        ncandidates = ncandidates + 1
        soslf = soslyapforward(s, p, dynamicfor(s, dyn))
        μ = measurefor(μs, dyn)
        cur = dot(μ, soslf)
        if cur > best
            best = cur
            best_dyn = dyn
        end
    end

    if ncandidates == 0
        error("$curstate does not have any incoming path of length $l.")
    end
    if best_dyn == nothing
        error("Oops, this should not happend, please report this bug.")
    end
    best_dyn
end

function sosbuilditeration(s::AbstractDiscreteSwitchedSystem, seq, μs, p_k, l, Δt, curstate, iter)
    best_dyn = best_dynamic(s, μs, p_k, l, curstate)

    curstate = state(s, best_dyn.seq[1], false)
    append!(seq, best_dyn)
    x = vars(p_k)
    p_k = p_k(best_dyn.A * x, x)
    iter+l, curstate, p_k(best_dyn.A * x, x)
end

function sosbuilditeration(s::AbstractContinuousSwitchedSystem, seq, μs, p_k, l, Δt, curstate, iter)
    best_dyn = best_dynamic(s, μs, p_k, l, curstate)

    if iter == 1 || seq.seq[iter-1][1] == best_dyn
        curstate = state(s, best_dyn, false)
        x = vars(p_k)
        iter+1, curstate, p_k(integratorfor(seq.seq[iter-1]) * x, x)
    else
        mode, dt = seq.seq[iter-1]
        seq.seq[iter-1] = (mode, dt/2)
        iter, curstate, p_k
    end
end

# Extracting trajectory from Lyapunov
function sosbuildsequence(s::DiscreteSwitchedSystem, d::Integer; v_0=:Random, p_0=:Random, l::Integer=1, Δt=.01, niter::Integer=42, solver=MathProgBase.defaultSDPsolver, tol=1e-5)
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
        curstate = rand(1:1)
    else
        if !isa(v_0, Integer) || v_0 < 1 || v_0 > 1
            throw(ArgumentError("Invalid v_0=$v_0"))
        end
        curstate = v_0
    end

    p_k = p_0
    n = dim(s)
    seq = SwitchingSequence(s, niter)

    iter = 1
    while iter <= niter
        iter, curstate, p_k = sosbuilditeration(s, seq, lyap.dual, p_k, l, Δt, curstate, iter)
    end
    @assert seq.len == length(seq.seq)

    smp = Nullable{DiscretePeriodicSwitching}()
    for i = 1:seq.len
        P = speye(n)
        startNode = state(s, seq.seq[i], false)
        k = 0
        for j = i:seq.len
            mode = seq.seq[j]
            k = k + nlabels(s, mode)
            P = integratorfor(s, mode) * P
            if state(s, mode, true) == startNode
                lambda = ρ(P)
                growthrate = adaptgrowthrate(abs(lambda), @view seq.seq[i:j])
                if isnull(smp) || isbetter(growthrate, length(i:j), get(smp))
                    smp = Nullable{DiscretePeriodicSwitching}(DiscretePeriodicSwitching(s, seq.seq[i:j], growthrate))
                end
            end
        end
    end

    if !isnull(smp)
        updatesmp!(s, get(smp))
    end
    smp
end

#function sosbuildsequence(s::DiscreteSwitchedSystem, d::Integer; solver=MathProgBase.defaultSDPsolver, tol=1e-5)
#end

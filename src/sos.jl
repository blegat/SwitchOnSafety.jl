using JuMP
using PolyJuMP
using SumOfSquares
using MathProgBase
# getdual of JuMP and MathProgBase.SolverInterface conflict
using MathProgBase.SolverInterface.AbstractMathProgSolver

export getlyap, soslyap, soslyapb, sosbuildsequence

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
function getlyap(s::AbstractSwitchedSystem, d::Int; solver::AbstractMathProgSolver=JuMP.UnsetSolver(), tol=1e-5)
    if d > length(s.lyaps) || isnull(s.lyaps[d])
        soslyapb(s, d, solver=solver, tol=1e-5, cached=true)
    end
    get(s.lyaps[d])
end

# Building the Lyapunov constraints
function soslyapforward(s::AbstractDiscreteSwitchedSystem, p::Polynomial, path)
    xin = vars(s, startnode(path))
    xout = vars(s, endnode(path))
    p(dynamicfor(s, path) * xin, xout)
end
function soslyapforward(s::AbstractContinuousSwitchedSystem, p::Polynomial, mode::Int)
    x = vars(p)
    dot(differentiate(p, x), dynamicfor(s, mode) * x)
end
soslyapscaling(s::AbstractDiscreteSwitchedSystem, γ, d) = γ^(2*d)
soslyapscaling(s::AbstractContinuousSwitchedSystem, γ, d) = 2*d*γ
function soslyapconstraint(s::AbstractSwitchedSystem, model::JuMP.Model, p, edge, d, γ)
    getid(x) = x.id
    @polyconstraint model soslyapforward(s, lyapforout(p, edge), edge) <= soslyapscaling(s, γ, d) * lyapforin(p, edge)
end
function soslyapconstraints(s::AbstractSwitchedSystem, model::JuMP.Model, p, d, γ)
    [soslyapconstraint(s, model, p, mode, d, γ) for mode in modes(s, 1)]
end
function soslyapconstraints(s::ConstrainedDiscreteSwitchedSystem, model::JuMP.Model, p, d, γ)
    [soslyapconstraint(s, model, p, edge, d, γ) for edge in edges(s.G)]
end

function buildlyap(model::JuMP.Model, x::Vector{PolyVar{true}}, d::Int)
    Z = monomials(x, 2*d)
    @polyvariable model p Z
    @polyconstraint model p >= sum(x.^(2*d))
    p
end
lyapforin(p::Vector, mode::Int) = p[1]
lyapforout(p::Vector, mode::Int) = p[1]
lyapforin(p::Vector, edge::Edge) = p[edge.first]
lyapforout(p::Vector, edge::Edge) = p[edge.second]

# Solving the Lyapunov problem
function soslyap(s::AbstractSwitchedSystem, d, γ; solver::AbstractMathProgSolver=JuMP.UnsetSolver())
    model = JuMP.Model(solver=solver)
    p = [buildlyap(model, vars(s, v), d) for v in 1:nnodes(s)]
    cons = soslyapconstraints(s, model, p, d, γ)
    # I suppress the warning "Not solved to optimality, status: Infeasible"
    status = solve(model, suppress_warnings=true)
    if status == :Optimal
        status, getvalue.(p), nothing
    elseif status == :Infeasible
        status, nothing, getdual.(cons)
    else
        status, nothing, nothing
    end
end

function getsoslyapinitub(s::AbstractDiscreteSwitchedSystem, d::Integer)
    #_, sosub = pradiusb(s, 2*d)
    #sosub
    Inf
end
function getsoslyapinitub(s::AbstractContinuousSwitchedSystem, d::Integer)
    Inf
end

function increaselb(s::AbstractDiscreteSwitchedSystem, lb, step)
    lb *= step
end

soschecktol(soslb, sosub) = sosub - soslb
soschecktol(s::AbstractDiscreteSwitchedSystem, soslb, sosub) = soschecktol(log(soslb), log(sosub))
soschecktol(s::AbstractContinuousSwitchedSystem, soslb, sosub) = soschecktol(soslb, sosub)

sosshift(s::AbstractDiscreteSwitchedSystem, b, shift) = exp(log(b) + shift)
sosshift(s::AbstractContinuousSwitchedSystem, b, shift) = b + shift

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
    η = min(ρA(s), binomial(n+d-1, d))
    sosub / η^(1/(2*d))
end
soslb2lb(s::AbstractContinuousSwitchedSystem, soslb, d) = -Inf

# Obtaining bounds with Lyapunov
function soslyapb(s::AbstractSwitchedSystem, d::Integer; solver::AbstractMathProgSolver=JuMP.UnsetSolver(), tol=1e-5, step=1, cached=true)
    # The SOS ub is greater than the JSR hence also greater than any of its lower bound
    soslb = s.lb
    sosub = getsoslyapinitub(s, d)
    primal = dual = nothing
    while soschecktol(s, soslb, sosub) > tol
        mid = sosmid(s, soslb, sosub, step)
        status, curprimal, curdual = soslyap(s, d, mid, solver=solver)
        if !(status in [:Optimal, :Unbounded, :Infeasible])
            midlb = sosmid(s, soslb, mid, step)
            statuslb, curprimallb, curduallb = soslyap(s, d, midlb, solver=solver)
            if statuslb in [:Optimal, :Unbounded, :Infeasible]
                mid = midlb
                status = statuslb
                curprimal = curprimallb
                curdual = curduallb
            else
                midub = sosmid(s, mid, sosub, step)
                statusub, curprimalub, curdualub = soslyap(s, d, midub, solver=solver)
                if statusub in [:Optimal, :Unbounded, :Infeasible]
                    mid = midub
                    status = statusub
                    curprimal = curprimalub
                    curdual = curdualub
                end
            end
        end
        if status == :Optimal || status == :Unbounded # FIXME Unbounded is for a Mosek bug
            if !(curprimal === nothing) # FIXME remove
                primal = curprimal
            end
            sosub = mid
        elseif status == :Infeasible
            dual = curdual
            soslb = mid
        else
            warn("Solver returned with status : $statuslb for γ=$midlb, $status for γ=$mid and $statusub for γ=$midub. Stopping bisection with $(soschecktol(s, soslb, sosub)) > $tol (= tol)")
            break
        end
    end
    if cached
        if primal === nothing
            status, primal, _ = soslyap(s, d, sosub, solver=solver)
            @assert status == :Optimal
            @assert !(primal === nothing)
        end
        if dual === nothing
            status, _, dual = soslyap(s, d, soslb, solver=solver)
            if status != :Infeasible
                soslb = sosshift(s, sosub, -0.99 * tol)
                status, _, dual = soslyap(s, d, soslb, solver=solver)
                @assert status == :Infeasible
                @assert !(dual === nothing)
            end
        end
        setlyap!(s, Lyapunov(d, soslb, dual, sosub, primal))
    end
    ub = sosub
    lb = soslb2lb(s, soslb, d)
    updateb!(s, lb, ub)
end

#function sosbuildsequence(s::DiscreteSwitchedSystem, d::Integer; solver=MathProgBase.defaultSDPsolver, tol=1e-5)
#end

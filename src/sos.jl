using JuMP
using SumOfSquares

export getlyap, soslyap, soslyapb, soslyapbs

# Storing the Lyapunov
function setlyap!(s, lyap::Lyapunov)
    d = lyap.d
    lyaps = getlyaps(s)
    if length(lyaps) < d
        sizehint!(lyaps, d)
        while length(lyaps) < d
            push!(lyaps, nothing)
        end
    end
    lyaps[d] = lyap
end
function getlyap(s::AbstractSwitchedSystem, d::Int; kws...)
    lyaps = getlyaps(s)
    if d > length(lyaps) || lyaps[d] === nothing
        soslyapb(s, d; cached=true, kws...)
    end
    return lyaps[d]
end

function getsoslyapinitub(s::AbstractDiscreteSwitchedSystem, d::Integer)
    #_, sosub = pradiusb(s, 2*d)
    #sosub
    return Inf
end
#function getsoslyapinitub(s::AbstractContinuousSwitchedSystem, d::Integer)
#    return Inf
#end

function getsoslyapinit(s, d)
    lyaps = getlyaps(s)
    if d <= length(lyaps) && lyaps[d] !== nothing
        lyap = lyaps[d]
        return lyap.soslb, lyap.dual, lyap.sosub, lyap.primal
    else
        # The SOS ub is greater than the JSR hence also greater than any of its lower bound.
        # Hence getlb(s) can be used as an initial lowerbound
        return getlb(s), nothing, getsoslyapinitub(s, d), nothing
    end
end

function superset(s::HybridSystems.AbstractHybridSystem, mode, d)
    q = GramMatrix(SOSDecomposition(variables(s, mode).^d))
    return SetProg.Sets.PolynomialSublevelSetAtOrigin(2d, q)
end

# For values of γ far from 1.0, it is better to divide A_i's by γ,
# it results in a problem that is better conditioned.
# This is clearly visible in [Example 5.4, PJ08] for which the JSR is ≈ 8.9
# So we ask the user to give a scaled system, i.e. do `p(A/γx) ≥ p(x)` rather
# than do `p(Ax) ≥ γ p(x)`
"""
    soslyap(s::AbstractSwitchedSystem, d; factory=nothing)

Find Sum-of-Squares Lyapunov functions; i.e. solves [(5), PJ08]
or gives moment matrices certifying the infeasibility of the problem.
Use [`ScaledHybridSystem`](@ref) to use a different growth rate than 1.

[PJ08] P. Parrilo and A. Jadbabaie.
*Approximation of the joint spectral radius using sum of squares*.
Linear Algebra and its Applications, Elsevier, **2008**, 428, 2385-2402
"""
function soslyap(s::HybridSystems.AbstractHybridSystem, d; factory=nothing)
    sets = HybridSystems.state_property(s, PolynomialLyapunov{Float64})
    infeasibility_certificates = HybridSystems.transition_property(s, MeasureLyapunov{Float64})
    set_variables = PolySet[PolySet(symmetric=true, degree=2d, superset=superset(s, mode, d))
                            for mode in states(s)]
    status = invariant_sets!(sets, 1:nstates(s), s, factory, set_variables,
                             volume_heuristic = nothing,
                             infeasibility_certificates = infeasibility_certificates,
                             verbose=0)
    if isinfeasible(status, true)
        @assert !isfeasible(status, true)
        if status[3] == MOI.NO_SOLUTION
            return status, nothing, nothing
        else
            return status, nothing, infeasibility_certificates
        end
    elseif isfeasible(status, true)
        status, sets, nothing
    else
        @assert !isdecided(status, true)
        status, nothing, nothing
    end
end

function increaselb(s::AbstractDiscreteSwitchedSystem, lb, step)
    lb *= step
end

soschecktol(soslb, sosub) = sosub - soslb
soschecktol(s::AbstractDiscreteSwitchedSystem, soslb, sosub) = soschecktol(log(soslb), log(sosub))
tol_diff_str(s::AbstractDiscreteSwitchedSystem) = "Log-diff   "
#soschecktol(s::AbstractContinuousSwitchedSystem, soslb, sosub) = soschecktol(soslb, sosub)

sosshift(s::AbstractDiscreteSwitchedSystem, b, shift) = exp(log(b) + shift)
function sosshift(s::AbstractDiscreteSwitchedSystem, b, shift, scaling)
    return sosshift(s, b / scaling, shift) * scaling
end
#sosshift(s::AbstractContinuousSwitchedSystem, b, shift) = b + shift

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
usestep(soslb, sosub) = isfinite(soslb) ⊻ isfinite(sosub)
sosmid(s::AbstractDiscreteSwitchedSystem, soslb, sosub, step) = exp(sosmid(log(soslb), log(sosub), step))
usestep(s::AbstractDiscreteSwitchedSystem, soslb, sosub) = usestep(log(soslb), log(sosub))
#sosmid(s::AbstractContinuousSwitchedSystem, soslb, sosub, step) = sosmid(soslb, sosub, step)
function sosmid(s::AbstractDiscreteSwitchedSystem, soslb, sosub, step, scaling)
    sosmid(s, soslb / scaling, sosub / scaling, step) * scaling
end

function soslb2lb(s::AbstractDiscreteSwitchedSystem, soslb, d)
    n = maximum(statedim.(s, states(s)))
    η = min(ρA(s), binomial(n+d-1, d))
    soslb / η^(1/(2*d))
end
#soslb2lb(s::AbstractContinuousSwitchedSystem, soslb, d) = -Inf

function showbs(s, soslb, sosub, tol, verbose, ok::Bool)
    if verbose >= 2 || (ok && verbose >= 1)
        println("Lower bound: $soslb")
        println("Upper bound: $sosub")
        println("$(tol_diff_str(s)): $(soschecktol(s, soslb, sosub)) $(ok ? '≤' : '>') $tol")
    end
end

function showmid(γ, status, verbose)
    if verbose >= 3
        println("  Trial value of γ: $γ")
        println("Termination status: $(status[1])")
        println("     Primal status: $(status[2])")
        println("       Dual status: $(status[3])")
        println("        Raw status: $(status[4])")
        if !isdecided(status, true)
            problem_status = "Unknown"
        elseif isfeasible(status, true)
            problem_status = "Feasible"
        else
            @assert isinfeasible(status, true)
            problem_status = "Infeasible"
        end
        println("    Problem status: $problem_status")
    end
end

# Binary Search

"""
    soslyapbs(s::AbstractSwitchedSystem, d::Integer,
              soslb, dual,
              sosub, primal;
              verbose=0, tol=1e-5, step=0.5, scaling=quickub(s),
              ranktols=tol, disttols=tol, kws...)


Find the smallest `γ` such that [`soslyap`](@ref) is feasible.
"""
function soslyapbs(s::AbstractSwitchedSystem, d::Integer,
                   soslb, dual,
                   sosub, primal;
                   verbose=0, tol=1e-5, step=0.5, scaling=quickub(s),
                   ranktols=tol, disttols=tol, kws...)
    _lyap(γ) = soslyap(ScaledHybridSystem(s, γ), d; kws...)
    while soschecktol(s, soslb, sosub) > tol
        showbs(s, soslb, sosub, tol, verbose, false)
        mid = sosmid(s, soslb, sosub, step, scaling)
        status, curprimal, curdual = _lyap(mid)
        showmid(mid, status, verbose)
        if !isdecided(status, true)
            if usestep(s, soslb, sosub)
                step *= 2
                continue
            end
            # If mid-tol/2 and mid+tol/2 also Stall, there would be an interval of length tol of Stall -> impossible to satisfy requirements
            # the distance between soslb and mid is at least tol/2.
            # Sometimes, mid is far from soslb and is at a point where the solver Stall even if it is far from the optimum point.
            # In that case, it is better to take (mid + soslb)/2
            midlb = min(sosmid(s, soslb, mid, step, scaling),
                        sosshift(s, mid, -tol/2, scaling))
            # If mid-tol/2 is too close to soslb, we would not make progress!
            # So we ensure we make a progress of at least tol/8. If dual is nothing, then that would still be progress to find a dual
            if dual !== nothing
                midlb = max(midlb, sosshift(s, soslb, tol/8))
            end
            statuslb, curprimallb, curduallb = _lyap(midlb)
            showmid(midlb, statuslb, verbose)
            if isdecided(statuslb, true)
                mid = midlb
                status = statuslb
                curprimal = curprimallb
                curdual = curduallb
            else
                midub = max(sosmid(s, mid, sosub, step, scaling),
                            sosshift(s, mid, tol/2, scaling))
                if primal !== nothing
                    midub = min(midub, sosshift(s, sosub, -tol/8))
                end
                statusub, curprimalub, curdualub = _lyap(midub)
                showmid(midub, statusub, verbose)
                if isdecided(statusub, true)
                    mid = midub
                    status = statusub
                    curprimal = curprimalub
                    curdual = curdualub
                end
            end
        end
        if isinfeasible(status, true)
            dual = curdual
            if dual !== nothing
                sosextractcycle(s, dual, d, ranktols=ranktols, disttols=disttols)
            end
            soslb = mid
        elseif isfeasible(status, true)
            if !(curprimal === nothing) # FIXME remove
                primal = curprimal
            end
            sosub = mid
        else
            @warn("Solver returned with status : $statuslb for γ=$midlb, $status for γ=$mid and $statusub for γ=$midub. Stopping bisection with $(soschecktol(s, soslb, sosub)) > $tol (= tol)")
            break
        end
    end
    if soschecktol(s, soslb, sosub) ≤ tol # it is not guaranteed because of the break
        showbs(s, soslb, sosub, tol, verbose, true)
    end
    soslb, dual, sosub, primal
end

"""
    soslyapbs(s::AbstractSwitchedSystem, d::Integer,
              soslb, dual,
              sosub, primal;
              verbose=0, tol=1e-5, step=0.5, scaling=quickub(s),
              ranktols=tol, disttols=tol, kws...)


Find upper bounds to the (constrained) Joint Spectral Radius [(5), PJ08].
Lower bounds a obtained using guarantees in [LPJ19].

[LPJ19] B. Legat, P. A. Parrilo and R. M. Jungers
*An entropy-based bound for the computational complexity of a switched system*.
IEEE Transactions on Automatic Control, **2019**.

[PJ08] P. Parrilo and A. Jadbabaie.
*Approximation of the joint spectral radius using sum of squares*.
Linear Algebra and its Applications, Elsevier, **2008**, 428, 2385-2402
"""
function soslyapb(s::AbstractSwitchedSystem, d::Integer; factory=nothing, tol=1e-5, cached=true, kws...)
    soslb, dual, sosub, primal = soslyapbs(s::AbstractSwitchedSystem, d::Integer, getsoslyapinit(s, d)...; factory=factory, tol=tol, kws...)
    _lyap(γ) = soslyap(ScaledHybridSystem(s, γ), d, factory=factory)
    if cached
        if primal === nothing
            if isfinite(sosub)
                status, primal, _ = _lyap(sosub)
                @assert isfeasible(status, true)
                @assert primal !== nothing
            else
                error("Bisection ended with infinite sosub=$sosub")
            end
        end
        if dual === nothing
            if isfinite(soslb)
                status, _, dual = _lyap(soslb)
                if !isinfeasible(status, true)
                    soslb = sosshift(s, soslb, -tol)
                    status, _, dual = _lyap(soslb)
                    if !isinfeasible(status, true)
                        @warn("We ignore getlb and start from scratch. tol was probably set too small and soslb is too close to the JSR so soslb-tol is too close to the JSR")
                        soslb = 0. # FIXME fix for continuous
                    end
                    soslb, dual, sosub, primal = soslyapbs(s::AbstractSwitchedSystem, d::Integer, soslb, dual, sosub, primal; factory=factory, tol=tol, kws...)
                    @assert dual !== nothing
                end
            else
                error("Bisection ended with infinite soslb=$soslb")
            end
        end
        setlyap!(s, Lyapunov(d, soslb, dual, sosub, primal))
    end
    ub = sosub
    lb = soslb2lb(s, soslb, d)
    updateb!(s, lb, ub)
end

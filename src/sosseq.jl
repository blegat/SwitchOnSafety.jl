export sosbuildsequence

function candidates(s::AbstractDiscreteSwitchedSystem, l, curstate)
    switchings(s, l, curstate, false)
end

function best_dynamic(s::AbstractSwitchedSystem, d, μs, p::GramMatrix, l, curstate, candidates)
    function rating(dyn, transformation)
        # The major part of the time spent in sosbuildsequence is done
        # in this line so the transformation is cached in `transformation`
        soslf = SetProg.apply_transformation(p, transformation)
        ν = measurefor(μs, dyn)
        # return dot(μ, soslf)
        return dot(getmat(ν), getmat(soslf))
    end
    best_dyn = nothing
    best = -Inf
    for (dyn, transformation) in candidates[curstate]
        cur = rating(dyn, transformation)
        if cur > best
            best = cur
            best_dyn = dyn
        end
    end
    if best_dyn === nothing
        error("$curstate does not have any incoming path of length $l.")
    end
    return best_dyn
end

function apply_map(sw::HybridSystems.DiscreteSwitchingSequence, p, new_vars, d)
    return SetProg.apply_matrix(p, sw.A, new_vars, d)
end

function sosbuilditeration(s::AbstractDiscreteSwitchedSystem, d, seq, μs, p_k::GramMatrix, l, Δt, curstate, iter, candidates)
    best_dyn = best_dynamic(s, d, μs, p_k, l, curstate, candidates)

    curstate = state(s, best_dyn.seq[1], false)
    prepend!(seq, best_dyn)
    x = variables(p_k)
    iter+l, curstate, apply_map(best_dyn, p_k, variables(s, source(s, best_dyn)), d)
end

#function sosbuilditeration(s::AbstractContinuousSwitchedSystem, seq, μs, p_prev, l, Δt, curstate, iter)
#    dyn = best_dynamic(s, μs, p_prev, l, curstate)
#    ub = Δt
#    x = variables(p_prev)
#    p_cur = p_prev(integratorfor(s, (dyn, ub)) * x, x)
#    best_dyn = best_dynamic(s, μs, p_cur, l, curstate)
#    if best_dyn != dyn
#        lb = 0
#        while ub-lb > 1e-5
#            mid = (lb+ub)/2
#            p_cur = p_prev(integratorfor(s, (dyn, mid)) * x, x)
#            best_dyn = best_dynamic(s, μs, p_cur, l, curstate)
#            if best_dyn == dyn
#                lb = mid
#            else
#                ub = mid
#            end
#        end
#    end
#
#    p_cur = p_prev(integratorfor(s, (dyn, ub)) * x, x)
#    if iter > 1 && dyn == seq.seq[iter-1][1] && duration(@view seq.seq[1:(iter-1)]) < 1000# && ub < 1e-3
#        seq.seq[iter-1] = (dyn, seq.seq[iter-1][2]+ub)
#    else
#        push!(seq, (dyn, ub))
#        curstate = state(s, dyn, false)
#        iter += 1
#    end
#    iter, curstate, p_cur
#end

function transformation(s::AbstractSwitchedSystem, path, d)
    xin = variables(s, source(s, path))
    xout = variables(s, target(s, path))
    return SetProg.transformation(monomials(xout, d), dynamicfort(s, path), xin, d)
end

# Extracting trajectory from Lyapunov

"""
    sosbuildsequence(s::AbstractSwitchedSystem, d::Integer;
                     v_0=:Random, p_0=:Random, l::Integer=1,
                     Δt::Float64=1., niter::Integer=42,
                     kws...)

Compute the truncation of length `l` of the high growth rate infinite sequence
produced by the algorithm introduced in [LJP17]. The trajectory ends at mode
`v_0` (or a random one if `v_0` is `:Random`) and is built backward as explained
in [LJP17]. The measures used to guide the construction are the infeasibility
certificates of highest growth rate computed by [`soslyap`](@ref) with
polynomials of degree `2d`. The starting polynomial is either `p_0`, or a random
strictly sum-of-squares polynomial if `p_0` is `:Random` or the primal solution
of [`soslyap`](@ref) certifying the best upper bound for mode `v_0`.

* [LJP17] B. Legat, R. M. Jungers, and P. A. Parrilo.
[Certifying unstability of Switched Systems using Sum of Squares Programming](https://arxiv.org/abs/1710.01814),
arXiv preprint arXiv:1710.01814, **2017**.
"""
function sosbuildsequence(s::AbstractSwitchedSystem, d::Integer;
                          v_0=:Random, p_0=:Random, l::Integer=1,
                          Δt::Float64=1., niter::Integer=42,
                          kws...)
    lyap = getlyap(s, d; kws...)

    if v_0 == :Random
        curstate = rand(states(s))
    else
        if !(v_0 in states(s))
            throw(ArgumentError("Invalid v_0=$v_0"))
        end
        curstate = v_0
    end

    if p_0 == :Primal
        p_k = lyap.primal[curstate].p
    elseif p_0 == :Random
        Z = monomials(variables(s, curstate), d)
        p_k = randsos(Z, monotype=:Gram, r=1)
    else
        # otherwise p_0 is assumed to be an sos polynomial given by the user
        p_k = p_0::GramMatrix
    end

    seq = HybridSystems.switchingsequence(s, niter, curstate)

    candidates = Dict(state => [(dyn, transformation(s, dyn, d)) for dyn in collect(switchings(s, l, state, false))]
                      for state in states(s))

    iter = 1
    while iter <= niter
        iter, curstate, p_k = sosbuilditeration(s, d, seq, lyap.dual, p_k, l, Δt, curstate, iter, candidates)
        # Avoid having it go to zero
        p_k = gram_operate(/, p_k, p_k(variables(p_k) => ones(Int, nvariables(p_k))))
    end
    @assert seq.len == length(seq.seq)
    seq
end

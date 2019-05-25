export invariant_polytopes
using LinearAlgebra

struct SymPoint{VT}
    point::VT
end

include("conitope.jl")
include("balanced_real_polytope.jl")
include("balanced_complex_polytope.jl")

const PolytopeLike = Union{Conitope, BalancedRealPolytope, BalancedComplexPolytope}

function _reset_model(p::PolytopeLike)
    p.model = nothing
    p.z = nothing
    p.t_0 = nothing
end
function Polyhedra.convexhull!(p::PolytopeLike, v::SymPoint)
    push!(p.points, v.point)
    # Invalidates the model
    _reset_model(p)
end
function _fix(p::PolytopeLike, v)
    JuMP.fix.(p.z, v)
end
function Base.in(v, brp::PolytopeLike)
    if isempty(brp.points)
        return false
    end
    if brp.model === nothing
        _parametrized_model(brp, v)
    else
        _fix(brp, v)
    end
    JuMP.optimize!(brp.model)
    status = termination_status(brp.model)
    if status == MOI.OPTIMAL
        return true
    elseif status == MOI.INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED
        return false
    else
        error("Solver returned status $status.")
    end
end
function in_ratio(v, p::PolytopeLike)
    if isempty(p.points)
        # The polytope is the set containing the origin only so `v` needs to be
        # multiplied by 0.0 to be mapped to the polytope
        return 0.0
    end
    _ratio_model(p, v)
    @objective(p.model, Max, p.t_0)
    JuMP.optimize!(p.model)
    status = termination_status(p.model)
    if status == MOI.OPTIMAL
        r = JuMP.value(p.t_0)
        _reset_model(p)
        return r
    else
        _reset_model(p)
        error("Solver returned status $status.")
    end
end

function root_path(parent, leaf, n)
    path = Int[]
    cur = leaf
    while cur > n
        push!(path, cur)
        cur = parent[cur]
    end
    return cur, reverse(path)
end

"""
Computes invariant balanced polytopes for system `s` using the algorithm of
[GP13] with the spectral maximizing product candidate (s.m.p.) `smp`. The
choice between Complex Balanced Polytope, Real Balanced Polytope or Conitope
using invariance of the positive orthant is done automatically. The solver used
is `factory`. This solver should support second order cone if the Complex
Balanced Polytope is used (e.g. if the leading eigenvalue has an imaginary
part). The maximum depth considered is `max_length`, a point `x` is considered
to be in the interior of the polytope if the maximum `λ` such that `λ * x` is
in the polytope is larger than `1 + tol`.

New s.m.p. candidates are considered with tolerance `new_candidate_tol`. A
maximum of `max_cycles` times the current s.m.p. candidate is prepended to this
new candidate to try to find a candidate with growth rate higher than the
current one, we also stop adding cycles if the length reaches `max_smp_length`.
If `gready`, then cycles continue to be prepended even if the it has a higher growth
rate than the current s.m.p. until the growth rate decreases when prepending a new
cycle.

If `verbose` is 1, the number of leaves is printed when the depth is a multiple
of `log_step_length`.

[GP13] N. Guglielmi and V. Protasov.
*Exact computation of joint spectral characteristics of linear operators*.
Foundations of Computational Mathematics 13.1, **2013**, 37-97.
"""
function invariant_polytopes(
    s::AbstractDiscreteSwitchedSystem, factory::JuMP.OptimizerFactory,
    smp::AbstractPeriodicSwitching; max_length=10, verbose=2, tol=nothing,
    max_cycles=10, new_candidate_tol=1e-6, gready=false, max_smp_length=50,
    log_step_length=10)
    leaves = Int[]
    sets = PolytopeLike[]
    new_smp = true
    max_lowtol = -Inf
    min_uptol = Inf
    max_low_candidate_tol = -Inf
    min_up_candidate_tol = Inf
    while new_smp
        new_smp = false
        n = length(smp.period)
        # The ith vertex is between transition[i] and transition[i+1]
        transition = [smp.period[end]; smp.period[1:(end-1)]]
        parent = [n; 1:(n-1)]
        A = prod(reverse(SwitchOnSafety.integratorfor.(s, smp.period)))
        # non-Hermitian matrices are not diagonalizable by StaticArrays
        if A isa SMatrix && !ishermitian(A)
            A = Matrix(A)
        end
        eig = eigen(Matrix(A))
        # TODO check whether the multiplicity is one
        λind = argmax(abs.(eig.values))
        λ_raw = eig.values[λind]
        λ = abs(λ_raw)^(1/n)
        @assert λ ≈ smp.growthrate
        v1 = eig.vectors[:, λind]
        u1 = inv(eig.vectors)[λind, :]
        if isreal(λ_raw)
            v1 = real.(v1)
            u1 = real.(u1)
        end
        if all(all(x -> x ≥ 0, SwitchOnSafety.integratorfor(s, t))
               for t in transitions(s))
            if any(x -> x < 0, v1)
                v1 = -v1
            end
            @assert all(x -> x ≥ 0, v1)
            sets = [Conitope(statedim(s, mode), Vector{eltype(v1)}[], factory)
                    for mode in modes(s)]
        else
            if isreal(λ_raw)
                sets = [BalancedRealPolytope(statedim(s, mode),
                                             Vector{eltype(v1)}[],
                                             factory)
                        for mode in modes(s)]
            else
                sets = [BalancedComplexPolytope(statedim(s, mode),
                                                Vector{eltype(v1)}[],
                                                factory)
                        for mode in modes(s)]
            end
        end
        vertices = [v1]
        duals = [u1 / dot(u1, v1)]
        for i in 2:n
            # TODO should maybe also add complex conjugate if it is complex like in CJSR
            push!(vertices,
                  SwitchOnSafety.integratorfor(s, transition[i]) * vertices[i-1] / λ)
        end
        for i in n:-1:1
            j = (i % n) + 1
            u = SwitchOnSafety.integratorfor(s, transition[j])' * vertices[j]
            u /= dot(u, vertices[i])
            push!(duals, u)
        end
        _mode(i) = target(s, transition[i])
        for i in 1:n
            if verbose ≥ 3
                println("v_$i = ", vertices[i])
            end
            Polyhedra.convexhull!(sets[_mode(i)], SymPoint(vertices[i]))
        end
        leaves = collect(1:n)
        for k in 1:max_length
            new_smp && break
            isempty(leaves) && break
            if verbose ≥ 2 || (verbose ≥ 1 && ((k % log_step_length) == 0))
                println("Depth $k: $(map(p -> length(p.points), sets)) points, $(length(leaves)) living leaves...")
            end
            new_leaves = Int[]
            for leaf in leaves
                new_smp && break
                root, path = root_path(parent, leaf, n)
                if verbose ≥ 3
                    println("Path starting starting at root $root: $(map(i -> symbol(s, transition[i]), path)):")
                end
                for t in out_transitions(s, _mode(leaf))
                    new_smp && break
                    if k == 1 && t == transition[(leaf % n) + 1]
                        continue
                    end
                    v = SwitchOnSafety.integratorfor(s, t) * vertices[leaf] / λ
                    if verbose ≥ 3
                        print(symbol(s, t), " : ", v, " : ")
                    end
                    if tol === nothing
                        is_in = v in sets[target(s, t)]
                    else
                        r = in_ratio(v, sets[target(s, t)]) - 1
                        is_in = r ≥ tol
                    end
                    if is_in
                        if verbose ≥ 3
                            printstyled("dead leaf", color=:red)
                            if tol !== nothing
                                println(": $r ≥ $tol")
                            end
                            println()
                        end
                        if tol !== nothing
                            min_uptol = min(min_uptol, r)
                        end
                    else
                        if target(s, t) == target(s, transition[root])
                            candidate_rating = abs(real(dot(duals[root], v))) - 1
                            if candidate_rating > new_candidate_tol
                                min_up_candidate_tol = min(min_up_candidate_tol, candidate_rating)
                                period = [[transition[i] for i in path]; t]
                                if verbose ≥ 2
                                    printstyled("New candidate s.m.p. found", bold=true, color=:blue)
                                    println(" because ⟨v_$root*, v⟩ - 1 = $(dot(duals[root], v) - 1) > $new_candidate_tol:")
                                end
                                cycle = [transition[root+1:n]; transition[1:root]]
                                _smp = periodicswitching(s, period, scaling = smp.growthrate)
                                __smp = _smp
                                better = isbetter(_smp, smp)
                                n_cycles = 0
                                for i in 1:max_cycles
                                    if length(_smp.period) + length(cycle) > max_smp_length
                                        break
                                    end
                                    if better
                                        if gready
                                            if i > 1 && !isbetter(_smp, __smp)
                                                # It does not improve, let's stop
                                                _smp = __smp
                                                break
                                            end
                                        else
                                            break
                                        end
                                    end
                                    if verbose ≥ 3
                                        print(_smp)
                                        if better
                                            print(" is better than existing s.m.p. but it might still be improved")
                                        else
                                            print(" is not better yet")
                                        end
                                        println(", adding a cycle.")
                                    end
                                    period = [cycle; period]
                                    __smp = _smp
                                    _smp = periodicswitching(s, period, scaling = smp.growthrate)
                                    n_cycles += 1
                                    if !better
                                        better = isbetter(_smp, smp)
                                    end
                                end
                                if better
                                    new_smp = true
                                    smp = _smp
                                    if verbose ≥ 1
                                        println(smp)
                                    end
                                else
                                    if verbose ≥ 2
                                        printstyled("aborting", bold=true, color=:red)
                                        print(" with smp candidate $_smp of length $(length(_smp.period)) after prefixing the suffix with $n_cycles times the current s.m.p.")
                                        @assert n_cycles ≤ max_cycles
                                        if n_cycles == max_cycles
                                            println(" Increase the value of the `max_cycles` keyword argument to go further.")
                                        else
                                            println(" Increase the value of the `max_smp_length` keyword argument to go further.")
                                        end
                                    end
                                end
                            else
                                max_low_candidate_tol = max(max_low_candidate_tol, candidate_rating)
                            end
                        end
                        if !new_smp
                            if verbose ≥ 3
                                printstyled("living leaf", color=:green)
                                if tol !== nothing
                                    println(": $r < $tol")
                                end
                                println()
                            end
                            if tol !== nothing
                                max_lowtol = max(max_lowtol, r)
                            end
                            push!(transition, t)
                            push!(parent, leaf)
                            push!(vertices, v)
                            push!(new_leaves, length(parent))
                            Polyhedra.convexhull!(sets[_mode(length(parent))], SymPoint(v))
                        end
                    end
                end
            end
            leaves = new_leaves
        end
    end
    if verbose ≥ 1
        if isempty(leaves)
            print("0 living leaves, $smp is a ")
            printstyled("spectral maximizing product", bold=true)
            println(" (s.m.p.)")
        else
            print("Still $(length(leaves)) living leaves, ")
            printstyled("aborting", bold=true, color=:red)
            println(" at depth $max_length. Increase the `max_length` keyword argument to go further.")
        end
        if tol === nothing
            println("Use the `tol` keyword argument to change the threshold to evaluate dead and living leaves.")
        else
            println("Use `tol` < $max_lowtol to eliminate at least one more living leaf.")
            println("Use `tol` > $min_uptol to keep at least one more dead leaf.")
        end
        if isfinite(max_low_candidate_tol)
            println("Use `new_candidate_tol` < $max_low_candidate_tol to consider at least one more candidate as s.m.p.")
        end
        if isfinite(min_up_candidate_tol)
            println("Use `new_candidate_tol` > $min_up_candidate_tol to eliminate at least one s.m.p. candidate.")
        end
    end
    return smp, isempty(leaves), sets
end

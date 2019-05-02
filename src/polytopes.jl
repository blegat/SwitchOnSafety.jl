export invariant_polytopes
using LinearAlgebra

struct SymPoint{VT}
    point::VT
end

include("conitope.jl")
include("balanced_real_polytope.jl")

const PolytopeLike = Union{Conitope, BalancedRealPolytope}

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
function Base.in(v, brp::PolytopeLike)
    if brp.model === nothing
        _parametrized_model(brp, v)
    else
        JuMP.fix.(brp.z, v)
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

# Algorithm of [GP13]
#
# [GP13] N. Guglielmi and V. Protasov.
# *Exact computation of joint spectral characteristics of linear operators*.
# Foundations of Computational Mathematics 13.1, **2013**, 37-97.
function invariant_polytopes(
    s::AbstractDiscreteSwitchedSystem, factory::JuMP.OptimizerFactory,
    smp::AbstractPeriodicSwitching; max_length=10, verbose=1, tol=nothing,
    max_cycles=10, new_candidate_tol=1e-6, gready=false)
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
                error("Balanced Complex Polytope not supported yet.")
            end
        end
        vertices = [v1]
        duals = [u1 / dot(u1, v1)]
        for i in 2:n
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
            if verbose ≥ 2
                println("v_$i = ", vertices[i])
            end
            Polyhedra.convexhull!(sets[_mode(i)], SymPoint(vertices[i]))
        end
        leaves = collect(1:n)
        for k in 1:max_length
            new_smp && break
            isempty(leaves) && break
            if verbose ≥ 1
                println("Depth $k: $(map(p -> length(p.points), sets)) points, $(length(leaves)) living leaves...")
            end
            new_leaves = Int[]
            for leaf in leaves
                new_smp && break
                root, path = root_path(parent, leaf, n)
                if verbose ≥ 2
                    println("Path starting starting at root $root: $(map(i -> symbol(s, transition[i]), path)):")
                end
                for t in out_transitions(s, _mode(leaf))
                    new_smp && break
                    if k == 1 && t == transition[(leaf % n) + 1]
                        continue
                    end
                    v = SwitchOnSafety.integratorfor(s, t) * vertices[leaf] / λ
                    if verbose ≥ 2
                        print(symbol(s, t), " : ", v, " : ")
                    end
                    if tol === nothing
                        is_in = v in sets[target(s, t)]
                    else
                        r = in_ratio(v, sets[target(s, t)]) - 1
                        is_in = r ≥ tol
                    end
                    if is_in
                        if verbose ≥ 2
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
                        candidate_rating = dot(duals[root], v) - 1
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
                            for i in 1:max_cycles
                                if better
                                    if gready
                                        if i > 1 && !isbetter(_smp, __smp)
                                            # It does not improve, let's stop
                                            break
                                        end
                                    else
                                        break
                                    end
                                end
                                if verbose ≥ 2
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
                                if !better
                                    better = isbetter(_smp, smp)
                                end
                            end
                            if better
                                new_smp = true
                                smp = _smp
                                if verbose ≥ 2
                                    println(smp)
                                end
                            else
                                if verbose ≥ 1
                                    printstyled("aborting", bold=true, color=:red)
                                    println(" with smp candidate $_smp after prefixing suffix with $max_cycles times the current s.m.p. Increase the `max_cycles` keyword argument to go further.")
                                end
                            end
                        else
                            max_low_candidate_tol = max(max_low_candidate_tol, candidate_rating)
                        end
                        if !new_smp
                            if verbose ≥ 2
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
            println("Use `tol` < $max_lowtol to eliminite at least one more living leaf.")
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

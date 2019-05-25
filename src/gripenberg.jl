export gripenberg

struct BranchState{MT<:AbstractMatrix, ET}
    A::MT
    seq::Vector{ET}
    mode::Int
    p::Float64 # min_j ||prod_{i=1}^j A_i||^{1/j}
end

"""
    gripenberg(s::AbstractDiscreteSwitchedSystem; δ=1e-2,
                    max_eval = 10000, max_ρ_eval = max_eval,
                    max_norm_eval = max_eval, max_length = 50,
                    matrix_norm = A -> opnorm(A, 2), verbose = 1)

Gripenberg algorithm [G96] for computing an upper bound `ub` and a lower bound
`lb` to the Joint Spectral Radius such that `ub - lb ≤ 1e-2`.

[G96] Gripenberg, G. Computing the joint spectral radius.
*Linear Algebra and its Applications*, *Elsevier*, **1996**, *234*, 43-60
"""
function gripenberg(s::AbstractDiscreteSwitchedSystem; δ=1e-2,
                    max_eval = 10000, max_ρ_eval = max_eval,
                    max_norm_eval = max_eval, max_length = 50,
                    matrix_norm = A -> opnorm(A, 2), verbose = 1)
    MT = typeof(dynamicfort(s, first(out_transitions(s, 1))))
    ET = transitiontype(s)
    branches = [BranchState{MT, ET}(Matrix(LinearAlgebra.I, HybridSystems.statedim(s, mode), HybridSystems.statedim(s, mode)), ET[], mode, Inf) for mode in modes(s)]
    smp = nothing
    ub = Inf
    n_ρ_eval = 0
    n_norm_eval = 0
    cur_length = 0
    while !isempty(branches) && n_ρ_eval < max_ρ_eval &&
        n_norm_eval < max_norm_eval && cur_length < max_length &&
        (smp === nothing || ub > smp.growthrate + δ)
        @assert cur_length == length(first(branches).seq)
        cur_length += 1
        new_branches = eltype(branches)[]
        for branch in branches
            for t in out_transitions(s, branch.mode)
                A = dynamicfort(s, t) * branch.A
                seq = [branch.seq; t]
                if source(s, first(seq)) == target(s, last(seq))
                    new_smp = periodicswitching(s, seq, A)
                    n_ρ_eval += 1
                    notifyperiodic!(s, new_smp)
                    if smp === nothing || isbetter(new_smp, smp)
                        smp = new_smp
                    end
                end
                p = min(matrix_norm(A)^(1/length(seq)), branch.p)
                n_norm_eval += 1
                if smp === nothing || p > smp.growthrate + δ
                    push!(new_branches, BranchState(A, seq, target(s, t), p))
                end
            end
        end
        branches = new_branches
        β = 0.0
        if !isempty(branches)
            β = max(β, maximum(branch -> branch.p, branches))
        end
        if smp !== nothing
            β = max(smp.growthrate + δ, β)
        end
        ub = min(ub, β)
    end
    if verbose ≥ 1
        function _keyword_feedback(cur, max)
            if cur < max
                print(cur, " < ", max)
            else
                printstyled(cur, " ≥ ", max, bold=true, color=:red)
            end
        end
        print("ρ evaluations   : ")
        _keyword_feedback(n_ρ_eval, max_ρ_eval)
        println(" = max_ρ_eval.")
        print("norm evaluations: ")
        _keyword_feedback(n_norm_eval, max_norm_eval)
        println(" = max_norm_eval.")
        print("switch length   : ")
        _keyword_feedback(cur_length, max_length)
        println(" = max_length.")
    end
    return smp, ub
end

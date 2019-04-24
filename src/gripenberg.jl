export gripenberg

struct BranchState{MT<:AbstractMatrix, ET}
    A::MT
    seq::Vector{ET}
    mode::Int
    p::Float64 # min_j ||prod_{i=1}^j A_i||^{1/j}
end

"""
   gripenberg(s; δ=1e-2, max_eval=1000, max_length = 50, matrix_norm=A -> opnorm(A, 2))

[G96] Gripenberg, G. Computing the joint spectral radius.
*Linear Algebra and its Applications*, *Elsevier*, **1996**, *234*, 43-60
"""
function gripenberg(s; δ=1e-2, max_eval=1000, max_length = 50, matrix_norm=A -> opnorm(A, 2))
    MT = typeof(SwitchOnSafety.dynamicfort(s, first(SwitchOnSafety.io_transitions(s, 1, true))))
    ET = transitiontype(s)
    branches = [BranchState{MT, ET}(Matrix(LinearAlgebra.I, HybridSystems.statedim(s, mode), HybridSystems.statedim(s, mode)), ET[], mode, Inf) for mode in modes(s)]
    smp = nothing
    ub = Inf
    n_eval = 0
    while !isempty(branches) && n_eval < max_eval && length(first(branches).seq) < max_length &&
        (smp === nothing || ub > smp.growthrate + δ)
        new_branches = eltype(branches)[]
        for branch in branches
            for t in SwitchOnSafety.io_transitions(s, mode, true)
                A = SwitchOnSafety.dynamicfort(s, t) * branch.A
                seq = [branch.seq; t]
                if source(s, first(seq)) == target(s, last(seq))
                    new_smp = SwitchOnSafety.periodicswitching(s, seq, A)
                    SwitchOnSafety.notifyperiodic!(s, new_smp)
                    if smp === nothing || SwitchOnSafety.isbetter(new_smp, smp)
                        smp = new_smp
                    end
                    n_eval += 1
                end
                p = min(matrix_norm(A)^(1/length(seq)), branch.p)
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
    return smp, ub
end

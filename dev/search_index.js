var documenterSearchIndex = {"docs":
[{"location":"invariant/#Invariant-Sets-1","page":"Invariant Sets","title":"Invariant Sets","text":"","category":"section"},{"location":"invariant/#","page":"Invariant Sets","title":"Invariant Sets","text":"Given a system described with MathematicalSystems.jl and HybridSystems.jl, an invariant set can be computed using the following function:","category":"page"},{"location":"invariant/#","page":"Invariant Sets","title":"Invariant Sets","text":"invariant_sets\ninvariant_sets!","category":"page"},{"location":"invariant/#SwitchOnSafety.invariant_sets","page":"Invariant Sets","title":"SwitchOnSafety.invariant_sets","text":"invariant_sets(system::AbstractHybridSystem, optimizer_constructor,\n               set_variables::AbstractVector{<:SetProg.AbstractVariable};\n               volume_heuristic = nth_root,\n               infeasibility_certificates = nothing,\n               verbose=1,\n               λ=Dict{transitiontype(system), Float64}(),\n               enabled = states(system))\n\nCompute maximal invariant sets of the family set_variables for the modes of system using the solver provided by optimizer_constructor. The volume of the sets is estimated using volume_heuristic. If the program is infeasible, the certificates for each transition are stored in infeasibility_certificates. For the containment of non-homogeneous, the S-procedure might be a Bilinear Matrix Inequality (BMI) which is NP-hard to solve. To avoid that, provides the value of λ to use in the dictionary λ. To ignore some state and the transitions involving these states in the computation, give an enabled vector without them.\n\n\n\n\n\n","category":"function"},{"location":"invariant/#SwitchOnSafety.invariant_sets!","page":"Invariant Sets","title":"SwitchOnSafety.invariant_sets!","text":"invariant_sets!(sets, modes_to_compute, system::AbstractHybridSystem,\n                args...; kwargs...)\n\nSimilar to invariant_sets(system, args...; kwargs...) but stores the result in sets and only compute the modes in modes_to_compute. The other sets in enabled that are not in modes_to_compute are assumed to have the value given in sets.\n\n\n\n\n\n","category":"function"},{"location":"#Switch-On-Safety-(SOS)-1","page":"Index","title":"Switch On Safety (SOS)","text":"","category":"section"},{"location":"#","page":"Index","title":"Index","text":"This packages implements methods for computing invariant sets using Sum Of Squares Programming. It supports:","category":"page"},{"location":"#","page":"Index","title":"Index","text":"Systems defined in MathematicalSystems.jl.\nHybrid Systems defined in HybridSystems.jl.","category":"page"},{"location":"#","page":"Index","title":"Index","text":"It also includes utilities for approximation the Joint Spectral Radius.","category":"page"},{"location":"#Contents-1","page":"Index","title":"Contents","text":"","category":"section"},{"location":"#","page":"Index","title":"Index","text":"Pages = [\"invariant.md\", \"jsr.md\"]\nDepth = 2","category":"page"},{"location":"jsr/#Joint-Spectral-Radius-(JSR)-1","page":"Joint Spectral Radius","title":"Joint Spectral Radius (JSR)","text":"","category":"section"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"The Joint Spectral Radius (JSR) of a discrete-time switched system is the minimal value of γ such that the system scaled by γ is asymptotically stable.","category":"page"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"ScaledHybridSystem","category":"page"},{"location":"jsr/#SwitchOnSafety.ScaledHybridSystem","page":"Joint Spectral Radius","title":"SwitchOnSafety.ScaledHybridSystem","text":"struct ScaledHybridSystem{T, H <: Union{DiscreteSwitchedLinearSystem, ConstrainedDiscreteSwitchedLinearSystem}} <: HybridSystems.AbstractHybridSystem\n    system::H\n    γ::T\nend\n\nDiscrete-time system where each reset map is scaled by γ, that is, the reset map x ↦ Ax of system is replaced by x ↦ Ax/γ.\n\n\n\n\n\n","category":"type"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"The JSR is NP-hard to compute but several methods exist to approximate this quantity.","category":"page"},{"location":"jsr/#Gripenberg-algorithm,-a-Branch-and-Bound-approach-1","page":"Joint Spectral Radius","title":"Gripenberg algorithm, a Branch-and-Bound approach","text":"","category":"section"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"The following approachs searches through all the different switching sequences, pruning some using the running estimate of the lower bound.","category":"page"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"gripenberg","category":"page"},{"location":"jsr/#SwitchOnSafety.gripenberg","page":"Joint Spectral Radius","title":"SwitchOnSafety.gripenberg","text":"gripenberg(s::AbstractDiscreteSwitchedSystem; δ=1e-2,\n                max_eval = 10000, max_ρ_eval = max_eval,\n                max_norm_eval = max_eval, max_length = 50,\n                matrix_norm = A -> opnorm(A, 2), verbose = 1)\n\nGripenberg algorithm [G96] for computing an upper bound ub and a lower bound lb to the Joint Spectral Radius such that ub - lb ≤ 1e-2.\n\n[G96] Gripenberg, G. Computing the joint spectral radius. Linear Algebra and its Applications, Elsevier, 1996, 234, 43-60\n\n\n\n\n\n","category":"function"},{"location":"jsr/#Sum-of-Squares-approach-1","page":"Joint Spectral Radius","title":"Sum-of-Squares approach","text":"","category":"section"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"The following method computes upper bounds using Sum Of Squares Programming","category":"page"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"soslyap\nsoslyapbs","category":"page"},{"location":"jsr/#SwitchOnSafety.soslyap","page":"Joint Spectral Radius","title":"SwitchOnSafety.soslyap","text":"soslyap(s::AbstractSwitchedSystem, d; optimizer_constructor=nothing)\n\nFind Sum-of-Squares Lyapunov functions; i.e. solves [(5), PJ08] or gives moment matrices certifying the infeasibility of the problem. Use ScaledHybridSystem to use a different growth rate than 1.\n\n[PJ08] P. Parrilo and A. Jadbabaie. Approximation of the joint spectral radius using sum of squares. Linear Algebra and its Applications, Elsevier, 2008, 428, 2385-2402\n\n\n\n\n\n","category":"function"},{"location":"jsr/#SwitchOnSafety.soslyapbs","page":"Joint Spectral Radius","title":"SwitchOnSafety.soslyapbs","text":"soslyapbs(s::AbstractSwitchedSystem, d::Integer,\n          soslb, dual,\n          sosub, primal;\n          verbose=0, tol=1e-5, step=0.5, scaling=quickub(s),\n          ranktols=tol, disttols=tol, kws...)\n\nFind the smallest γ such that soslyap is feasible.\n\n\n\n\n\n","category":"function"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"The infeasibility certificates computed in the binary search carried out by soslyapbs can be used to produce cycles of high growth rate using findsmp on the sequence produced by sosbuildsequence.","category":"page"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"sosbuildsequence\nfindsmp","category":"page"},{"location":"jsr/#SwitchOnSafety.sosbuildsequence","page":"Joint Spectral Radius","title":"SwitchOnSafety.sosbuildsequence","text":"sosbuildsequence(s::AbstractSwitchedSystem, d::Integer;\n                 v_0=:Random, p_0=:Random, l::Integer=1,\n                 Δt::Float64=1., niter::Integer=42,\n                 kws...)\n\nCompute the truncation of length l of the high growth rate infinite sequence produced by the algorithm introduced in [LJP17]. The trajectory ends at mode v_0 (or a random one if v_0 is :Random) and is built backward as explained in [LJP17]. The measures used to guide the construction are the infeasibility certificates of highest growth rate computed by soslyap with polynomials of degree 2d. The starting polynomial is either p_0, or a random strictly sum-of-squares polynomial if p_0 is :Random or the primal solution of soslyap certifying the best upper bound for mode v_0.\n\n[LJP17] B. Legat, R. M. Jungers, and P. A. Parrilo.\n\nCertifying unstability of Switched Systems using Sum of Squares Programming, arXiv preprint arXiv:1710.01814, 2017.\n\n\n\n\n\n","category":"function"},{"location":"jsr/#SwitchOnSafety.findsmp","page":"Joint Spectral Radius","title":"SwitchOnSafety.findsmp","text":"findsmp(seq::HybridSystems.DiscreteSwitchingSequence)\n\nExtract the cycle of highest growth rate in the sequence seq.\n\n\n\n\n\n","category":"function"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"Alternatively, sosextractcycle can be used to find cycles of high growth rate.","category":"page"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"sosextractcycle","category":"page"},{"location":"jsr/#SwitchOnSafety.sosextractcycle","page":"Joint Spectral Radius","title":"SwitchOnSafety.sosextractcycle","text":"sosextractcycle(s::AbstractDiscreteSwitchedSystem, dual, d::Integer;\n                ranktols=1e-5, disttols=1e-5)\n\nExtract cycles of high growth rate from atomic occupation measures given by the infeasibility certificates of highest growth rate computed by soslyap. The method is detailed in [LJP17].\n\n[LJP17] B. Legat, R. M. Jungers, and P. A. Parrilo.\n\nCertifying unstability of Switched Systems using Sum of Squares Programming, arXiv preprint arXiv:1710.01814, 2017.\n\n\n\n\n\n","category":"function"},{"location":"jsr/#Polytopic-approach-1","page":"Joint Spectral Radius","title":"Polytopic approach","text":"","category":"section"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"The following method can verify numerically that a cycle is an s.m.p. by computing invariant polytopes. When it is not an s.m.p., it find a candidate of higher growth rate.","category":"page"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"@docs` invariant_polytopes","category":"page"}]
}

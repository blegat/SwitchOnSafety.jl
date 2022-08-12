var documenterSearchIndex = {"docs":
[{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"EditURL = \"https://github.com/blegat/SwitchOnSafety.jl/blob/master/examples/PJ08e54.jl\"","category":"page"},{"location":"generated/PJ08e54/#Example-5.4-of-[PJ08]-1","page":"PJ08e54","title":"Example 5.4 of [PJ08]","text":"","category":"section"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"(Image: ) (Image: ) Contributed by: Benoît Legat","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"Example 5.4 of [PJ08]. The JSR is conjectured to be equal to 89149 = sqrtrho(A_1 A_3) in [PJ08] by noticing that rho_textSOS-4 approx rho_textSOS-6 approx 892.","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"In this notebook, we see that the technique developped in [LPJ17] finds the cycle (1 3) for degrees 1, 2, 3 already with l = 1. The upper bound is computed with Mosek v8.1.0.61 with an log-accuracy of 4e-7 (i.e. the different of the logarithms is smaller than 4 times 10^-7). We see that rho_textSOS-2  sqrtrho(A_1 A_3) as already shown in [PJ08] but we show that rho_textSOS-4  sqrtrho(A_1 A_3) and that rho_textSOS-6 - sqrtrho(A_1 A_3)  153 times 10^-7 hence conjecture that they are equal.","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"[PJ08] P. Parrilo and A. Jadbabaie, Approximation of the joint spectral radius using sum of squares. Linear Algebra and its Applications, Elsevier, 2008, 428, 2385-2402","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"[LPJ17] Legat, B., Parrilo, P. A., & Jungers, R. M. Certifying unstability of Switched Systems using Sum of Squares Programming. arXiv preprint arXiv:1710.01814, 2017.","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"using SwitchOnSafety\nA1 = [ 0  1  7  4;\n       1  6 -2 -3;\n      -1 -1 -2 -6;\n       3  0  9  1]\nA2 = [-3  3  0 -2;\n      -2  1  4  9;\n       4 -3  1  1;\n       1 -5 -1 -2]\nA3 = [ 1  4  5 10;\n       0  5  1 -4;\n       0 -1  4  6;\n      -1  5  0  1]\ns = discreteswitchedsystem([A1, A2, A3])","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"We first apply Gripenberg","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"@time psw, ub = gripenberg(s)","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"Pick an SDP solver from this list.","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"import CSDP\noptimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);\n\nsosdata(s).lb = 0.0\n@time lb2, ub2 = soslyapb(s, 1, optimizer_constructor=optimizer_constructor, tol=4e-7, step=0.5, verbose=1)","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"From the infeasibility certificate of the last infeasible SemiDefinite Program (SDP) solved by in order to obtain the upper bound, we find the lower bound rho(A_1A_3)^12 approx 8914964.","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"seq = sosbuildsequence(s, 1)\npsw = findsmp(seq)","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"We see below that rho_textSOS-4  8919  sqrtrho(A_1 A_3).","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"sosdata(s).lb = 0.0\n@time lb4, ub4 = soslyapb(s, 2, optimizer_constructor=optimizer_constructor, tol=4e-7, step=0.5, verbose=1)","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"We now try to round the infeasibility certificate to get a lower bound.","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"seq = sosbuildsequence(s, 2)\npsw = findsmp(seq)","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"We now compute rho_textSOS-6.","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"sosdata(s).lb = 0.0\n@time lb6, ub6 = soslyapb(s, 3, optimizer_constructor=optimizer_constructor, tol=4e-7, step=0.5, verbose=1)","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"To see how good is the upper bound, we use rounding to get a lower bound.","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"seq = sosbuildsequence(s, 3)\npsw = findsmp(seq)","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"We see below that rho_textSOS-6 - sqrtrho(A_1 A_3)  153 times 10^-7","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"ub6 - psw.growthrate","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"","category":"page"},{"location":"generated/PJ08e54/#","page":"PJ08e54","title":"PJ08e54","text":"This page was generated using Literate.jl.","category":"page"},{"location":"invariant/#Invariant-Sets-1","page":"Invariant Sets","title":"Invariant Sets","text":"","category":"section"},{"location":"invariant/#","page":"Invariant Sets","title":"Invariant Sets","text":"Given a system described with MathematicalSystems.jl and HybridSystems.jl, an invariant set can be computed using the following function:","category":"page"},{"location":"invariant/#","page":"Invariant Sets","title":"Invariant Sets","text":"invariant_sets\ninvariant_sets!","category":"page"},{"location":"invariant/#SwitchOnSafety.invariant_sets","page":"Invariant Sets","title":"SwitchOnSafety.invariant_sets","text":"invariant_sets(system::AbstractHybridSystem, optimizer_constructor,\n               set_variables::AbstractVector{<:SetProg.AbstractVariable};\n               volume_heuristic = nth_root,\n               infeasibility_certificates = nothing,\n               verbose=1,\n               λ=Dict{transitiontype(system), Float64}(),\n               enabled = states(system))\n\nCompute maximal invariant sets of the family set_variables for the modes of system using the solver provided by optimizer_constructor. The volume of the sets is estimated using volume_heuristic. If the program is infeasible, the certificates for each transition are stored in infeasibility_certificates. For the containment of non-homogeneous, the S-procedure might be a Bilinear Matrix Inequality (BMI) which is NP-hard to solve. To avoid that, provides the value of λ to use in the dictionary λ. To ignore some state and the transitions involving these states in the computation, give an enabled vector without them.\n\n\n\n\n\n","category":"function"},{"location":"invariant/#SwitchOnSafety.invariant_sets!","page":"Invariant Sets","title":"SwitchOnSafety.invariant_sets!","text":"invariant_sets!(sets, modes_to_compute, system::AbstractHybridSystem,\n                args...; kwargs...)\n\nSimilar to invariant_sets(system, args...; kwargs...) but stores the result in sets and only compute the modes in modes_to_compute. The other sets in enabled that are not in modes_to_compute are assumed to have the value given in sets.\n\n\n\n\n\n","category":"function"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"EditURL = \"https://github.com/blegat/SwitchOnSafety.jl/blob/master/examples/AJPR14e54.jl\"","category":"page"},{"location":"generated/AJPR14e54/#Example-5.4-of-[AJPR14]-1","page":"AJPR14e54","title":"Example 5.4 of [AJPR14]","text":"","category":"section"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"(Image: ) (Image: ) Contributed by: Benoît Legat","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"In the Example 5.4 of [AJPR14], the authors study the following system whose JSR is 3.917384715148. In this example, we reproduce their numerical results. The system is characterized by the matrices mathcalA = A_1 A_2 where the matrices are given below:","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"[AJPR14] A. Ahmadi, R. Jungers, P. Parrilo and M. Roozbehani, Joint spectral radius and path-complete graph Lyapunov functions. SIAM J. CONTROL OPTIM 52(1), 687-717, 2014.","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"using SwitchOnSafety\nA1 = [-1 -1\n      -4  0]\nA2 = [ 3  3\n      -2  1]\ns = discreteswitchedsystem([A1, A2])","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"Pick an SDP solver from this list.","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"import CSDP\noptimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);\nnothing #hide","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"We first try with a Common Quadratic Lyapunov Function (CQLF) and obtain rho_mathcalV^2(mathcalA) approx 3980613 and the corresponding lower bound rho_mathcalV^2(mathcalA)  sqrt2  2814547.","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"lb, ub = soslyapb(s, 1, optimizer_constructor=optimizer_constructor, tol=1e-6, verbose=1)","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"From the infeasibility certificate of the last infeasible SemiDefinite Program (SDP) solved by in order to obtain the upper bound, we find the lower bound rho(A_2A_1)^12 approx 3917384715148.","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"seq = sosbuildsequence(s, 1, p_0=:Primal)\npsw = findsmp(seq)","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"Now with a quartic lyapunov function, we get the tighter upper bound rho_mathcalV^textSOS4(mathcalA) approx 39241.","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"lb, ub = soslyapb(s, 2, optimizer_constructor=optimizer_constructor, tol=1e-6, verbose=1)","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"We find the same lower bound from the infeasibility certificate.","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"seq = sosbuildsequence(s, 2, p_0=:Primal)\npsw = findsmp(seq)","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"This upper bound is further improved by computing quadratic lyapunov functions for the the system obtained with the graph G_1 of [Fig. 3.1, AJPR14].","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"G1 = mdependentlift(s, 1, false)\nlb, ub = soslyapb(G1, 1, optimizer_constructor=optimizer_constructor, tol=1e-6, verbose=1)","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"We again find the same lower bound from the infeasibility certificate.","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"seq = sosbuildsequence(s, 2, p_0=:Primal)\npsw = findsmp(seq)","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"","category":"page"},{"location":"generated/AJPR14e54/#","page":"AJPR14e54","title":"AJPR14e54","text":"This page was generated using Literate.jl.","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"EditURL = \"https://github.com/blegat/SwitchOnSafety.jl/blob/master/examples/AP12e21.jl\"","category":"page"},{"location":"generated/AP12e21/#Example-2.1-of-[AP12]-1","page":"AP12e21","title":"Example 2.1 of [AP12]","text":"","category":"section"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"(Image: ) (Image: ) Contributed by: Benoît Legat","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"The JSR is 1.","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"[AP12] A. Ahmadi, and P. Parrilo Joint spectral radius of rank one matrices and the maximum cycle mean problem. CDC, 731-733, 2012.","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"using SwitchOnSafety\nA1 = [0 1\n      0 0]\nA2 = [0 0\n      1 0]\ns = discreteswitchedsystem([A1, A2])","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"Pick an SDP solver from this list.","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"import CSDP\noptimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);\nnothing #hide","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"We can see that we obtain the value of the JSR as upper and lower bound already with a CQLF. This is to be expected since the problem can be reduced to a maximum cycle mean problem which is easy to solve [AP12].","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"lb, ub = soslyapb(s, 1, optimizer_constructor=optimizer_constructor, tol=1e-4)","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"From the infeasibility certificate of the last infeasible SemiDefinite Program (SDP) solved by in order to obtain the upper bound, we find the lower bound rho(A_1A_2)^12 approx 1.","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"seq = sosbuildsequence(s, 1, p_0=:Primal)\npsw = findsmp(seq)","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"In fact, this same maximum cycle mean problem can be used to find the lower bound certificate as the infeasibility certificate of gamma close to 1 returned by the SDP solver contains atomic measures. The discrete problem formulated on the atoms is a maximum cycle means problem which gives the smp.","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"psw = sosextractcycle(s, 1, ranktols=1e-4)","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"","category":"page"},{"location":"generated/AP12e21/#","page":"AP12e21","title":"AP12e21","text":"This page was generated using Literate.jl.","category":"page"},{"location":"#Switch-On-Safety-(SOS)-1","page":"Index","title":"Switch On Safety (SOS)","text":"","category":"section"},{"location":"#","page":"Index","title":"Index","text":"This packages implements methods for computing invariant sets using Sum Of Squares Programming. It supports:","category":"page"},{"location":"#","page":"Index","title":"Index","text":"Systems defined in MathematicalSystems.jl.\nHybrid Systems defined in HybridSystems.jl.","category":"page"},{"location":"#","page":"Index","title":"Index","text":"It also includes utilities for approximation the Joint Spectral Radius.","category":"page"},{"location":"#Contents-1","page":"Index","title":"Contents","text":"","category":"section"},{"location":"#","page":"Index","title":"Index","text":"Pages = [\"invariant.md\", \"jsr.md\"]\nDepth = 2","category":"page"},{"location":"jsr/#Joint-Spectral-Radius-(JSR)-1","page":"Joint Spectral Radius","title":"Joint Spectral Radius (JSR)","text":"","category":"section"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"The Joint Spectral Radius (JSR) of a discrete-time switched system is the minimal value of γ such that the system scaled by γ is asymptotically stable.","category":"page"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"ScaledHybridSystem","category":"page"},{"location":"jsr/#SwitchOnSafety.ScaledHybridSystem","page":"Joint Spectral Radius","title":"SwitchOnSafety.ScaledHybridSystem","text":"struct ScaledHybridSystem{T, H <: Union{DiscreteSwitchedLinearSystem, ConstrainedDiscreteSwitchedLinearSystem}} <: HybridSystems.AbstractHybridSystem\n    system::H\n    γ::T\nend\n\nDiscrete-time system where each reset map is scaled by γ, that is, the reset map x ↦ Ax of system is replaced by x ↦ Ax/γ.\n\n\n\n\n\n","category":"type"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"The JSR is NP-hard to compute but several methods exist to approximate this quantity.","category":"page"},{"location":"jsr/#Gripenberg-algorithm,-a-Branch-and-Bound-approach-1","page":"Joint Spectral Radius","title":"Gripenberg algorithm, a Branch-and-Bound approach","text":"","category":"section"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"The following approachs searches through all the different switching sequences, pruning some using the running estimate of the lower bound.","category":"page"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"gripenberg","category":"page"},{"location":"jsr/#SwitchOnSafety.gripenberg","page":"Joint Spectral Radius","title":"SwitchOnSafety.gripenberg","text":"gripenberg(s::AbstractDiscreteSwitchedSystem; δ=1e-2,\n                max_eval = 10000, max_ρ_eval = max_eval,\n                max_norm_eval = max_eval, max_length = 50,\n                matrix_norm = A -> opnorm(A, 2), verbose = 1)\n\nGripenberg algorithm [G96] for computing an upper bound ub and a lower bound lb to the Joint Spectral Radius such that ub - lb ≤ 1e-2.\n\n[G96] Gripenberg, G. Computing the joint spectral radius. Linear Algebra and its Applications, Elsevier, 1996, 234, 43-60\n\n\n\n\n\n","category":"function"},{"location":"jsr/#Sum-of-Squares-approach-1","page":"Joint Spectral Radius","title":"Sum-of-Squares approach","text":"","category":"section"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"The following method computes upper bounds using Sum Of Squares Programming","category":"page"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"soslyap\nsoslyapbs","category":"page"},{"location":"jsr/#SwitchOnSafety.soslyap","page":"Joint Spectral Radius","title":"SwitchOnSafety.soslyap","text":"soslyap(s::AbstractSwitchedSystem, d; optimizer_constructor=nothing)\n\nFind Sum-of-Squares Lyapunov functions; i.e. solves [(5), PJ08] or gives moment matrices certifying the infeasibility of the problem. Use ScaledHybridSystem to use a different growth rate than 1.\n\n[PJ08] P. Parrilo and A. Jadbabaie. Approximation of the joint spectral radius using sum of squares. Linear Algebra and its Applications, Elsevier, 2008, 428, 2385-2402\n\n\n\n\n\n","category":"function"},{"location":"jsr/#SwitchOnSafety.soslyapbs","page":"Joint Spectral Radius","title":"SwitchOnSafety.soslyapbs","text":"soslyapbs(s::AbstractSwitchedSystem, d::Integer,\n          soslb, dual,\n          sosub, primal;\n          verbose=0, tol=1e-5, step=0.5, scaling=quickub(s),\n          ranktols=tol, disttols=tol, kws...)\n\nFind the smallest γ such that soslyap is feasible.\n\n\n\n\n\n","category":"function"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"The infeasibility certificates computed in the binary search carried out by soslyapbs can be used to produce cycles of high growth rate using findsmp on the sequence produced by sosbuildsequence.","category":"page"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"sosbuildsequence\nfindsmp","category":"page"},{"location":"jsr/#SwitchOnSafety.sosbuildsequence","page":"Joint Spectral Radius","title":"SwitchOnSafety.sosbuildsequence","text":"sosbuildsequence(s::AbstractSwitchedSystem, d::Integer;\n                 v_0=:Random, p_0=:Random, l::Integer=1,\n                 Δt::Float64=1., niter::Integer=42,\n                 kws...)\n\nCompute the truncation of length l of the high growth rate infinite sequence produced by the algorithm introduced in [LJP17]. The trajectory ends at mode v_0 (or a random one if v_0 is :Random) and is built backward as explained in [LJP17]. The measures used to guide the construction are the infeasibility certificates of highest growth rate computed by soslyap with polynomials of degree 2d. The starting polynomial is either p_0, or a random strictly sum-of-squares polynomial if p_0 is :Random or the primal solution of soslyap certifying the best upper bound for mode v_0.\n\n[LJP17] B. Legat, R. M. Jungers, and P. A. Parrilo.\n\nCertifying unstability of Switched Systems using Sum of Squares Programming, arXiv preprint arXiv:1710.01814, 2017.\n\n\n\n\n\n","category":"function"},{"location":"jsr/#SwitchOnSafety.findsmp","page":"Joint Spectral Radius","title":"SwitchOnSafety.findsmp","text":"findsmp(seq::HybridSystems.DiscreteSwitchingSequence)\n\nExtract the cycle of highest growth rate in the sequence seq.\n\n\n\n\n\n","category":"function"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"Alternatively, sosextractcycle can be used to find cycles of high growth rate.","category":"page"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"sosextractcycle","category":"page"},{"location":"jsr/#SwitchOnSafety.sosextractcycle","page":"Joint Spectral Radius","title":"SwitchOnSafety.sosextractcycle","text":"sosextractcycle(s::AbstractDiscreteSwitchedSystem, dual, d::Integer;\n                ranktols=1e-5, disttols=1e-5)\n\nExtract cycles of high growth rate from atomic occupation measures given by the infeasibility certificates of highest growth rate computed by soslyap. The method is detailed in [LJP17].\n\n[LJP17] B. Legat, R. M. Jungers, and P. A. Parrilo.\n\nCertifying unstability of Switched Systems using Sum of Squares Programming, arXiv preprint arXiv:1710.01814, 2017.\n\n\n\n\n\n","category":"function"},{"location":"jsr/#Polytopic-approach-1","page":"Joint Spectral Radius","title":"Polytopic approach","text":"","category":"section"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"The following method can verify numerically that a cycle is an s.m.p. by computing invariant polytopes. When it is not an s.m.p., it find a candidate of higher growth rate.","category":"page"},{"location":"jsr/#","page":"Joint Spectral Radius","title":"Joint Spectral Radius","text":"@docs` invariant_polytopes","category":"page"}]
}

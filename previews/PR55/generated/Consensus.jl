using LinearAlgebra
P = [ 1.0  0.0  1.0  2.0
     -1.0  0.0  1.0  2.0
      0.0  1.0  0.0 -3.0
      0.0 -1.0  0.0 -3.0
      0.0  0.0 -2.0  2.0]
for i in 1:4
    P[:, i] /= norm(P[:, i])
end
@show P' * ones(5)
round.(P' * P, digits=16)

using SwitchOnSafety
H1 = [0.2 0.8 0   0   0
      0   0.2 0.8 0   0
      0   0   0.2 0.8 0
      0   0   0   0.2 0.8
      0.2 0   0   0   0.8]
A1 = P' * H1 * P
H2 = [1   0   0   0   0
      0.8 0.2 0   0   0
      0   0   1   0   0
      0   0   0.8 0.2 0
      0   0   0   0   1]
A2 = P' * H2 * P
H3 = [1   0   0   0   0
      0   1   0   0   0
      0   0.8 0.2 0   0
      0   0   0   1   0
      0   0   0   0.2 0.8]
A3 = P' * H3 * P
function automaton(N)
    a = GraphAutomaton(N) # See [P17, Figure 2.17] for what automaton(3) should be
    add_transition!(a, 1, 1, 1) # Node i means, H1 was used i-1 steps ago
    for i in 2:N
        add_transition!(a, i-1, i, 2)
        add_transition!(a, i-1, i, 3)
        add_transition!(a, i, 1, 1)
    end
    return a
end
ss(N) = discreteswitchedsystem([A1, A2, A3], automaton(N))

import CSDP
optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);

ρ(A1)

soslyapb(ss(1), 1; optimizer_constructor=optimizer_constructor, tol=1e-7, verbose=1)

ss2 = ss(2);

soslyapb(ss2, 1; optimizer_constructor=optimizer_constructor, tol=1e-7, verbose=1)

@time seq = sosbuildsequence(ss2, 1, niter=100, l=1, p_0=:Primal)
@time psw = findsmp(seq)

soslyapb(ss2, 2; optimizer_constructor=optimizer_constructor, tol=1e-5, verbose=1)

@time seq = sosbuildsequence(ss2, 2, niter=100, l=2, p_0=:Primal)
@time psw = findsmp(seq)

soslyapb(ss2, 3; optimizer_constructor=optimizer_constructor, tol=1e-3, verbose=1)

@time seq = sosbuildsequence(ss2, 3, niter=45, l=1, p_0=:Primal)
@time psw = findsmp(seq)

soslyapb(ss(3), 1; optimizer_constructor=optimizer_constructor, tol=1e-7, verbose=1)

soslyapb(ss(4), 1; optimizer_constructor=optimizer_constructor, tol=1e-7, verbose=1)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl


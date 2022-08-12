# # Consensus

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Consensus.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Consensus.ipynb)
# **Contributed by**: Benoît Legat

# In this notebook, we show how to apply the JSR theory to compute the rate of convergence of agents to consensus.
# This example is [Example 2.52, P17].
#
# [P17] M. Philippe.
# *Path-Complete Methods and Analysis of Constrained Switching Systems*
# Doctoral dissertation, UCLouvain, **2017**

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

# We build the switched system as follows:

using Test #src
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
    # See [P17, Figure 2.17] for what automaton(3) should be
    a = GraphAutomaton(N)
    # Node i means, H1 was used i-1 steps ago
    add_transition!(a, 1, 1, 1)
    for i in 2:N
        add_transition!(a, i-1, i, 2)
        add_transition!(a, i-1, i, 3)
        add_transition!(a, i, 1, 1)
    end
    return a
end
ss(N) = discreteswitchedsystem([A1, A2, A3], automaton(N))

# Pick an SDP solver from [this list](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers).

import CSDP
optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);

# # $N = 1$
#
# With $N = 1$, the system is not switched, there is only one matrix and the JSR is approximately $0.7273$.

ρ(A1)

# With common quadratic Lyapunov functions, we get:

soslyapb(ss(1), 1; optimizer_constructor=optimizer_constructor, tol=1e-7, verbose=1)

# # $N = 2$

ss2 = ss(2);

# We start with CQLF (Common Quadratic Lyapunov Function)

soslyapb(ss2, 1; optimizer_constructor=optimizer_constructor, tol=1e-7, verbose=1)

# We use rounding to obtain good lower bounds:

@time seq = sosbuildsequence(ss2, 1, niter=100, l=1, p_0=:Primal)
@time psw = findsmp(seq)

# We now try with quartic polynomials

soslyapb(ss2, 2; optimizer_constructor=optimizer_constructor, tol=1e-5, verbose=1)

# Use rounding:

@time seq = sosbuildsequence(ss2, 2, niter=100, l=2, p_0=:Primal)
@time psw = findsmp(seq)

# Sextic:

soslyapb(ss2, 3; optimizer_constructor=optimizer_constructor, tol=1e-3, verbose=1)

# Rounding:

@time seq = sosbuildsequence(ss2, 3, niter=45, l=1, p_0=:Primal)
@time psw = findsmp(seq)

# # $N = 3$

soslyapb(ss(3), 1; optimizer_constructor=optimizer_constructor, tol=1e-7, verbose=1)

# # $N = 4$

soslyapb(ss(4), 1; optimizer_constructor=optimizer_constructor, tol=1e-7, verbose=1)

# # Example 5.4 of [AJPR14]

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/AJPR14e54.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/AJPR14e54.ipynb)
# **Contributed by**: Benoît Legat

# In the Example 5.4 of [AJPR14], the authors study the following system whose JSR is 3.917384715148.
# In this example, we reproduce their numerical results.
# The system is characterized by the matrices $\mathcal{A} = \{A_1, A_2\}$ where the matrices are given below:
#
# [AJPR14] A. Ahmadi, R. Jungers, P. Parrilo and M. Roozbehani,
# *Joint spectral radius and path-complete graph Lyapunov functions.*
# SIAM J. CONTROL OPTIM 52(1), 687-717, **2014**.

using Test #src
using SwitchOnSafety
A1 = [-1 -1
      -4  0]
A2 = [ 3  3
      -2  1]
s = discreteswitchedsystem([A1, A2])

# Pick an SDP solver from [this list](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers).

import CSDP
optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);

# We first try with a Common Quadratic Lyapunov Function (CQLF) and obtain $\rho_{\mathcal{V}^2}(\mathcal{A}) \approx 3.980613$ and the corresponding lower bound $\rho_{\mathcal{V}^2}(\mathcal{A}) / \sqrt{2} ≈ 2.814547$.

lb, ub = soslyapb(s, 1, optimizer_constructor=optimizer_constructor, tol=1e-6, verbose=1)
@test ub ≈ 3.9805 rtol=1e-4 #src

# From the infeasibility certificate of the last infeasible SemiDefinite Program (SDP) solved by in order to obtain the upper bound, we find the lower bound $\rho(A_2A_1)^{1/2} \approx 3.917384715148$.

seq = sosbuildsequence(s, 1, p_0=:Primal)
psw = findsmp(seq)
@test psw.growthrate ≈ 3.917384715148 #src

# Now with a quartic lyapunov function, we get the tighter upper bound
# $\rho_{\mathcal{V}^{\text{SOS},4}}(\mathcal{A}) \approx 3.9241$.

lb, ub = soslyapb(s, 2, optimizer_constructor=optimizer_constructor, tol=1e-6, verbose=1)
@test ub ≈ 3.9241 rtol=1e-4 #src

# We find the same lower bound from the infeasibility certificate.

seq = sosbuildsequence(s, 2, p_0=:Primal)
psw = findsmp(seq)
@test psw.growthrate ≈ 3.917384715148 #src

# This upper bound is further improved by computing quadratic lyapunov functions
# for the the system obtained with the graph $G_1$ of [Fig. 3.1, AJPR14].

G1 = mdependentlift(s, 1, false)
lb, ub = soslyapb(G1, 1, optimizer_constructor=optimizer_constructor, tol=1e-6, verbose=1)
@test ub ≈ 3.9224 rtol=1e-4 #src

# We again find the same lower bound from the infeasibility certificate.

seq = sosbuildsequence(s, 2, p_0=:Primal)
psw = findsmp(seq)
@test psw.growthrate ≈ 3.917384715148 #src

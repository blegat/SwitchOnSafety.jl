# # Example 2.1 of [AP12]

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/AJPR14e54.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/AJPR14e54.ipynb)
# **Contributed by**: Benoît Legat

# The JSR is 1.
#
# [AP12] A. Ahmadi, and P. Parrilo
# *Joint spectral radius of rank one matrices and the maximum cycle mean problem.*
# CDC, 731-733, **2012**.

using Test #src
using SwitchOnSafety
A1 = [0 1
      0 0]
A2 = [0 0
      1 0]
s = discreteswitchedsystem([A1, A2])

# Pick an SDP solver from [this list](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers).

import CSDP
optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);

# We can see that we obtain the value of the JSR as upper and lower bound already with a CQLF. This is to be expected since the problem can be reduced to a maximum cycle mean problem which is easy to solve [AP12].

lb, ub = soslyapb(s, 1, factory=factory, tol=1e-4)
@test ub ≈ 1 rtol=1e-4 #src

# From the infeasibility certificate of the last infeasible SemiDefinite Program (SDP) solved by in order to obtain the upper bound, we find the lower bound $\rho(A_2A_1)^{1/2} \approx 1$.

seq = sosbuildsequence(s, 1, p_0=:Primal)
psw = findsmp(seq)
@test psw.growthrate ≈ 1 #src

# In fact, this same maximum cycle mean problem can be used to find the lower bound certificate as the infeasibility certificate of $\gamma$ close to 1 returned by the SDP solver contains atomic measures. The discrete problem formulated on the atoms is a maximum cycle means problem which gives the smp.

psw = sosextractcycle(s, 1, ranktols=1e-4)
@test psw.growthrate ≈ 1 #src

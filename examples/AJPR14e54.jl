# # Example 5.4 of [AJPR14]

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/AJPR14e54.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/AJPR14e54.ipynb)
# **Contributed by**: Benoît Legat

# The JSR is 3.917384715148.
#
# [AJPR14] A. Ahmadi, R. Jungers, P. Parrilo and M. Roozbehani,
# *Joint spectral radius and path-complete graph Lyapunov functions.*
# SIAM J. CONTROL OPTIM 52(1), 687-717, **2014**.

using SwitchOnSafety
A1 = [-1 -1
      -4 0]
A2 = [ 3 3
      -2 1]
s = discreteswitchedsystem([A1, A2])

# Pick an SDP solver from [this list](http://jump.dev/JuMP.jl/v0.21.3/installation/#Getting-Solvers-1).

import CSDP
optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);

# We first try with CQLF...

soslyapb(s, 1, optimizer_constructor=optimizer_constructor, tol=1e-4)

# ... and try to find an smp from the infeasibility certificate of the SDP.

seq = sosbuildsequence(s, 1, p_0=:Primal)
psw = findsmp(seq)

# Now with quartic forms the upper bound get closer.

soslyapb(s, 2, optimizer_constructor=optimizer_constructor, tol=1e-4)

# The lower bound is the same.

seq = sosbuildsequence(s, 2, p_0=:Primal)
psw = findsmp(seq)

# # Example 5.4 of [PJ08]

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/PJ08e54.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/PJ08e54.ipynb)
# **Contributed by**: Benoît Legat

# Example 5.4 of [PJ08]. The JSR is conjectured to be equal to $8.9149 = \sqrt{\rho(A_1 A_3)}$ in [PJ08] by noticing that $\rho_{\text{SOS-}4} \approx \rho_{\text{SOS-}6} \approx 8.92$.
#
# In this notebook, we see that the technique developped in [LPJ17] finds the cycle $(1, 3)$ for degrees 1, 2, 3 already with $l = 1$.
# The upper bound is computed with Mosek v8.1.0.61 with an log-accuracy of `4e-7` (i.e. the different of the logarithms is smaller than $4 \times 10^{-7}$).
# We see that $\rho_{\text{SOS-}2} > \sqrt{\rho(A_1 A_3)}$ as already shown in [PJ08] but we show that $\rho_{\text{SOS-}4} > \sqrt{\rho(A_1 A_3)}$ and that $\rho_{\text{SOS-}6} - \sqrt{\rho(A_1 A_3)} < 1.53 \times 10^{-7}$ hence conjecture that they are equal.
#
# [PJ08] P. Parrilo and A. Jadbabaie,
# *Approximation of the joint spectral radius using sum of squares*.
# Linear Algebra and its Applications, Elsevier, **2008**, 428, 2385-2402
#
# [LPJ17] Legat, B., Parrilo, P. A., & Jungers, R. M.
# *Certifying unstability of Switched Systems using Sum of Squares Programming*.
# arXiv preprint arXiv:1710.01814, **2017**.

using Test #src
using SwitchOnSafety
A1 = [ 0  1  7  4;
       1  6 -2 -3;
      -1 -1 -2 -6;
       3  0  9  1]
A2 = [-3  3  0 -2;
      -2  1  4  9;
       4 -3  1  1;
       1 -5 -1 -2]
A3 = [ 1  4  5 10;
       0  5  1 -4;
       0 -1  4  6;
      -1  5  0  1]
s = discreteswitchedsystem([A1, A2, A3])

# We first apply Gripenberg

@time psw, ub = gripenberg(s)
@test ub ≈ 8.924964 rtol=1e-6 #src
@test psw.growthrate ≈ 8.914964 rtol=1e-6 #src

# Pick an SDP solver from [this list](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers).

import CSDP
optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);

sosdata(s).lb = 0.0
@time lb2, ub2 = soslyapb(s, 1, optimizer_constructor=optimizer_constructor, tol=4e-7, step=0.5, verbose=1)
@test ub2 ≈ 9.760675 rtol=1e-4 #src

# From the infeasibility certificate of the last infeasible SemiDefinite Program (SDP) solved by in order to obtain the upper bound, we find the lower bound $\rho(A_1A_3)^{1/2} \approx 8.914964$.

seq = sosbuildsequence(s, 1)
psw = findsmp(seq)
@test psw.growthrate ≈ 8.914964 rtol=1e-4 #src

# We see below that $\rho_{\text{SOS-}4} > 8.919 > \sqrt{\rho(A_1 A_3)}$.

sosdata(s).lb = 0.0
@time lb4, ub4 = soslyapb(s, 2, optimizer_constructor=optimizer_constructor, tol=4e-7, step=0.5, verbose=1)
@test ub4 ≈ 8.919820 rtol=1e-4 #src

# We now try to round the infeasibility certificate to get a lower bound.

seq = sosbuildsequence(s, 2)
psw = findsmp(seq)
@test psw.growthrate ≈ 8.914964 rtol=1e-4 #src

# We now compute $\rho_{\text{SOS-}6}$.

sosdata(s).lb = 0.0
@time lb6, ub6 = soslyapb(s, 3, optimizer_constructor=optimizer_constructor, tol=4e-7, step=0.5, verbose=1)
@test ub6 ≈ 8.914964 rtol=1e-4 #src

# To see how good is the upper bound, we use rounding to get a lower bound.

seq = sosbuildsequence(s, 3)
psw = findsmp(seq)
@test psw.growthrate ≈ 8.914964 rtol=1e-4 #src

# We see below that $\rho_{\text{SOS-}6} - \sqrt{\rho(A_1 A_3)} < 1.53 \times 10^{-7}$

ub6 - psw.growthrate

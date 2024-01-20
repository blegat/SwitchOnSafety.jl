using SwitchOnSafety
A1 = [-1 -1
      -4  0]
A2 = [ 3  3
      -2  1]
s = discreteswitchedsystem([A1, A2])

import CSDP
optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);

lb, ub = soslyapb(s, 1, optimizer_constructor=optimizer_constructor, tol=1e-6, verbose=1)

seq = sosbuildsequence(s, 1, p_0=:Primal)
psw = findsmp(seq)

lb, ub = soslyapb(s, 2, optimizer_constructor=optimizer_constructor, tol=1e-6, verbose=1)

seq = sosbuildsequence(s, 2, p_0=:Primal)
psw = findsmp(seq)

G1 = mdependentlift(s, 1, false)
lb, ub = soslyapb(G1, 1, optimizer_constructor=optimizer_constructor, tol=1e-6, verbose=1)

seq = sosbuildsequence(s, 2, p_0=:Primal)
psw = findsmp(seq)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

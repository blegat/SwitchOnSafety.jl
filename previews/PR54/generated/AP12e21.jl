using SwitchOnSafety
A1 = [0 1
      0 0]
A2 = [0 0
      1 0]
s = discreteswitchedsystem([A1, A2])

import CSDP
optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);

lb, ub = soslyapb(s, 1, optimizer_constructor=optimizer_constructor, tol=1e-4)

seq = sosbuildsequence(s, 1, p_0=:Primal)
psw = findsmp(seq)

psw = sosextractcycle(s, 1, ranktols=1e-4)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl


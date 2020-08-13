using SwitchOnSafety
A1 = [-1 -1
      -4 0]
A2 = [ 3 3
      -2 1]
s = discreteswitchedsystem([A1, A2])

import CSDP
optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);

soslyapb(s, 1, optimizer_constructor=optimizer_constructor, tol=1e-4)

seq = sosbuildsequence(s, 1, p_0=:Primal)
psw = findsmp(seq)

soslyapb(s, 2, optimizer_constructor=optimizer_constructor, tol=1e-4)

seq = sosbuildsequence(s, 2, p_0=:Primal)
psw = findsmp(seq)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl


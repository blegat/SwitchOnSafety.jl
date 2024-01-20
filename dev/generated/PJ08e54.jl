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

@time psw, ub = gripenberg(s)

import CSDP
optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true);

sosdata(s).lb = 0.0
@time lb2, ub2 = soslyapb(s, 1, optimizer_constructor=optimizer_constructor, tol=4e-7, step=0.5, verbose=1)

seq = sosbuildsequence(s, 1)
psw = findsmp(seq)

sosdata(s).lb = 0.0
@time lb4, ub4 = soslyapb(s, 2, optimizer_constructor=optimizer_constructor, tol=4e-7, step=0.5, verbose=1)

seq = sosbuildsequence(s, 2)
psw = findsmp(seq)

sosdata(s).lb = 0.0
@time lb6, ub6 = soslyapb(s, 3, optimizer_constructor=optimizer_constructor, tol=4e-7, step=0.5, verbose=1)

seq = sosbuildsequence(s, 3)
psw = findsmp(seq)

ub6 - psw.growthrate

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

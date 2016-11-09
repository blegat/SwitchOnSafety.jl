# The JSR is 3.917384715148
using SwitchedSystems
using Mosek
s = SwitchedSystem([[-1 -1; -4 0],[3 3; -2 1]])
#@show sosduallyapb(s, 1, solver=Mosek.MosekSolver(LOG=0))
@show sosduallyapb(s, 3, solver=Mosek.MosekSolver(LOG=0))

# [1] A. Ahmadi, R. Jungers, P. Parrilo and M. Roozbehani,
#     "Joint spectral radius and path-complete graph Lyapunov functions."
#     SIAM J. CONTROL OPTIM 52(1), 687-717, 2014.
#     The JSR is 3.917384715148

# Discrete Linear Unconstrained Switched Systems
dlu = [
       SwitchedSystem([[-1 -1; -4 0],[3 3; -2 1]]), # Example 5.4 of [1]

# Example 5.4 of
# A. Ahmadi, R. Jungers, P. Parrilo and M. Roozbehani,
# "Joint spectral radius and path-complete graph Lyapunov functions."
# SIAM J. CONTROL OPTIM 52(1), 687-717, 2014.
# The JSR is 3.917384715148

facts("Example 5.4") do
  ss = SwitchedSystem([[-1 -1; -4 0],[3 3; -2 1]])
  qub = 4.3195961
  @fact quicklb(ss) --> roughly(3)
  @fact quickub(ss) --> roughly(qub)
  qb = quickb(ss)
  @fact typeof(qb) --> NTuple{2,Float64}
  @fact qb[1] --> roughly(3)
  @fact qb[2] --> roughly(qub)
end

# Example 5.4 of
# A. Ahmadi, R. Jungers, P. Parrilo and M. Roozbehani,
# "Joint spectral radius and path-complete graph Lyapunov functions."
# SIAM J. CONTROL OPTIM 52(1), 687-717, 2014.
# The JSR is 3.917384715148

facts("Example 5.4") do
  s = SwitchedSystem([[-1 -1; -4 0],[3 3; -2 1]])
  qub = 4.3195961
  @fact quicklb(s) --> roughly(3)
  @fact quickub(s) --> roughly(qub)
  qb = quickb(s)
  @fact typeof(qb) --> NTuple{2,Float64}
  @fact qb[1] --> roughly(3)
  @fact qb[2] --> roughly(qub)
  # pnorm=Inf: [17.2409,7.42377,5.65853,4.92192,4.52597,4.28099,4.11275,3.99139,3.8991,3.827,3.76891,3.72121,3.6813,3.64744,3.61834,3.59307,3.57092,3.55135,3.53393,3.51832,3.50426,3.49152,3.47994]
  # pnorm=1:   [15.1245,7.35044,5.47402,4.84722,4.45847,4.23242,4.0709,3.95644,3.8687,3.80015,3.74485,3.69944,3.66141,3.62914,3.6014,3.5773,3.55616,3.53748,3.52085,3.50595,3.49253,3.48037,3.4693,3.45919]
  # pnorm=2:   [14.1535,6.81871,5.32651,4.7096,4.37067,4.15723,4.01117,3.90488,3.82407,3.7606,3.70945,3.66735,3.6321,3.60215,3.57639,3.554,3.53437,3.517,3.50154,3.48768,3.47519]
  @fact_throws ArgumentError pradius(s, 1)
  @fact_throws ArgumentError pradius(s, 2, :MagicAlgo)
  @fact pradius(s, 2) --> roughly(3.234535151)
  @fact pradius(s, 2, :BruteForce, pnorm=2) --> roughly(3.632097274822649)
  lb, ub = pradiusb(s, 2)
  @fact lb --> roughly(3.234535151)
  @fact ub --> roughly(4.574323479)
  lb, ub = pradiusb(s, 4)
  @fact lb --> roughly(3.427560156)
  @fact ub --> roughly(4.076078925)
  for solver in sdp_solvers
    context("With solver $(typeof(solver))") do
      lb, ub = sosb(s, 1, solver=solver)
      @fact lb --> roughly(2.814640557, rtol=1e-5)
      @fact ub --> roughly(3.980502849, rtol=1e-5)
      lb, ub = sosb(s, 2, solver=solver)
      @fact lb --> roughly(3.299750624, rtol=1e-5)
      @fact ub --> roughly(3.924086919, rtol=1e-5)
    end
  end
  @fact s.lb --> roughly(3.427560156)
  if isempty(sdp_solvers)
    @fact s.ub --> roughly(4.076078925)
  else
    @fact s.ub --> roughly(3.924086919, rtol=1e-5)
  end
end

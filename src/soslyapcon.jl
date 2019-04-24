    set_in = SetProg.Sets.PolynomialSublevelSetAtOrigin(d, p[source(s, edge)])
    set_out = SetProg.Sets.PolynomialSublevelSetAtOrigin(d, p[target(s, edge)])
    A = dynamicfort(s, edge, γ)
    @constraint(model, A * set_in ⊆ set_out)

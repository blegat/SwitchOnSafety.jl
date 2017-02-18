# Section 4 of
# M. Philippe, R. Essick, G. E. Dullerud and R. M. Jungers.
# "Stability of discrete-time switching systems with constrained switching sequences."
# Automatica, 72:242-250, 2016
# The JSR is 0.97482

const expected_lb = [0.7095,
                     0.7494,
                     0.8041,
                     0.8425,
                     0.8671,
                     0.8776]
const expected_ub = [1.0037,
                     0.9863,
                     0.9769,
                     0.9750,
                     0.9754,
                     0.9748]

@testset "[PEDJ] Section 4" begin
    A = [0.94 0.56; 0.14 0.46]
    B = [0; 1]
    k1 = -0.49
    k2 = 0.27
    As = [A + B * [k1 k2], A + B * [k1 0], A + B * [0 k2], A]
    G = DiGraph(4)
    σ = Dict{Edge, Int}()
    function add_edge_labeled!(G, u, v, σuv)
        e = u => v
        σ[e] = σuv
        add_edge!(G, e)
    end
    add_edge_labeled!(G, 1, 2, 3)
    add_edge_labeled!(G, 2, 1, 2)
    add_edge_labeled!(G, 1, 3, 1)
    add_edge_labeled!(G, 3, 1, 2)
    add_edge_labeled!(G, 2, 3, 1)
    add_edge_labeled!(G, 3, 2, 3)
    add_edge_labeled!(G, 3, 3, 1)
    add_edge_labeled!(G, 3, 4, 4)
    add_edge_labeled!(G, 4, 3, 1)
    s = ConstrainedDiscreteSwitchedSystem(As, G, σ)
    for solver in sdp_solvers
        println("  > With solver $(typeof(solver))")
        for d in 1:6
            tol = ismosek(solver) ? 1e-5 : 1e-4
            lb, ub = soslyapb(s, d, solver=solver, tol=tol)
            @show lb
            @show ub
        end
    end
end

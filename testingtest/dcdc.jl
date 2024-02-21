using Revise
using SwitchOnSafety

function controlswitch(A::AbstractVector{<:AbstractMatrix}, B::AbstractVector{<:AbstractMatrix})
    G = OneStateAutomaton(length(A))
    n = 1
    modes = ContinuousIdentitySystem.(map(s -> HybridSystems._getstatedim(A, G, s), states(G)))
    rm = LinearControlMap.(A, B)
    sw = HybridSystems.Fill(AutonomousSwitching(), n)
    HybridSystem(G, modes, rm, sw, Dict{Symbol, Any}())
end

r_l = 2
r_c = 2
r_0 = 2
x_l = 0.1
x_c = 3

a1_11 = (-(r_l / x_l))
a1_22 = (-((1 / x_c) * (1 / (r_0 + r_c))))
a2_11 = ((-(1 / x_l)) * (r_l + (((r_0 * r_c))/(r_0 + r_c))))
a2_12 = ((-(1 / x_l)) * (r_0 / (r_0 + r_c)))
a2_21 = ((1 / x_c) * (r_0 / (r_0 + r_c)))
a2_22 = ((-(1 / x_c)) * (1 / (r_0 + r_c)))

Δt = 1e-1

A1 = [
    a1_11 0 1/x_l
    0 a1_22 0
    0 0 0
]
A2 = [
    a2_11 a2_12 1 / x_l
    a2_21 a2_22 0
    0 0 0
]

B = reshape([0.0, 0.0, 1.0], 3, 1)

s = controlswitch([I + Δt * A1, I + Δt * A2], [B, B])
psw, ub = gripenberg(s, δ=0.004)

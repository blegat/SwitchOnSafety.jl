using Test
using SwitchOnSafety

function controlswitch(A::AbstractVector{<:AbstractMatrix}, B::AbstractVector{<:AbstractMatrix})
    G = OneStateAutomaton(length(A))
    n = 1
    modes = ContinuousIdentitySystem.(map(s -> HybridSystems._getstatedim(A, G, s), states(G)))
    rm = LinearControlMap.(A, B)
    sw = HybridSystems.Fill(AutonomousSwitching(), n)
    HybridSystem(G, modes, rm, sw, Dict{Symbol, Any}())
end

@testset "JCSR" begin
    A1 = [2 8; 2 5]
    B1 = [1 0; 0 1]
    A2 = [12 0; 10 0]
    B2 = [0 0; 1 0]
    s = controlswitch([A1, A2], [B1, B2])
    gripenberg(s)
end

r_l = 888
r_c = 203
r_0 = 333
x_l = 3
x_c = 5

a1_11 = (-(r_l / x_l))
a1_22 = (-((1 / x_c) * (1 / (r_0 + r_c))))
a2_11 = ((-(1 / x_l)) * (r_l + (((r_0 * r_c))/(r_0 + r_c))))
a2_12 = ((-(1 / x_l)) * (r_0 / (r_0 + r_c)))
a2_21 = ((1 / x_c) * (r_0 / (r_0 + r_c)))
a2_22 = ((-(1 / x_c)) * (1 / (r_0 + r_c)))


A1 = [1 0; 0 1] + [a1_11 0; 0 a1_22]
B1 = [1.0 0.0; 0.0 0.0]
A2 = [1 0; 0 1] + [a2_11 a2_12; a2_21 a2_22]
B2 = [1.0 0.0; 0.0 0.0]

s = controlswitch([A1, A2], [B1, B2])
controlswitch
gripenberg(s, δ=0.002)



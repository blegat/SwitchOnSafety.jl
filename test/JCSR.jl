using Revise
using SwitchOnSafety

m1 = LinearControlMap([2 2; 3 5], [1 0; 0 0])
m2 = LinearControlMap([1 2; 4 3], [0 1; 0 0])

function controlswitch(A::AbstractVector{<:AbstractMatrix}, B::AbstractVector{<:AbstractMatrix})
    G = OneStateAutomaton(length(A))
    n = 1
    modes = ContinuousIdentitySystem.(map(s -> HybridSystems._getstatedim(A, G, s), states(G)))
    rm = LinearControlMap.(A, B)
    sw = HybridSystems.Fill(AutonomousSwitching(), n)
    HybridSystem(G, modes, rm, sw, Dict{Symbol, Any}())
end

A1 = [2 8; 2 5]
B1 = [1 0; 0 1]
A2 = [12 0; 10 0]
B2 = [0 0; 1 0]
s = controlswitch([A1, A2], [B1, B2])
controlswitch
gripenberg(s)


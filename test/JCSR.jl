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


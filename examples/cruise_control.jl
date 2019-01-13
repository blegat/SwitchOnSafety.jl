using HybridSystems
include(joinpath(dirname(dirname(pathof(HybridSystems))), "examples", "cruise_control.jl"));

function constant(scaling)
    va = 15.6 * scaling
    vb = 24.5 * scaling
    vc = 29.5 * scaling
    v = (va, vb, vc)
    U = 4.0 * scaling
    D = 0.5 * scaling
    vmin = 5.0 * scaling
    return D, U, vmin, v
end
D, U, vmin, v = constant(1.0)
D = 1.0
U = 1.0
vmin = 1.0
v = (3.0, 4.0, 5.0)
T = 2
N = 1 + (T+1) * length(v)
m = 1.0
m0 = 1.0 #500;
h = 0.5 #0.8;
function system(M, T; v=v, N = 1 + (T+1) * length(v), sym=false, vmax = 35.)
    #return cruise_control_example(N, M, vmin = vmin, vmax=vmax, v=v, U=U, H=H, D=D, T=T, sym=sym, m0 = m0);
    return cruise_control_example(N, M, vmin = 1.0, vmax=vmax, v=v,
        U=1.0, H=h*T, D=1.0, T=T, sym=sym, m0 = 1.0, ks=1.0, m=1.0, kd=1.0)
end
_vec(M, d, v, u) = [repeat([d, v], M); v; u]
#Δv = (v[1] - vmin) / 2
Δv = (15.6 - 5) / 2
# |string_elongation| < D
# |acceleration| < U
sym_rect(M) = _vec(M, D, 1.0, U) # _vec(M, D, Δv, U)

using SetProg
using SwitchOnSafety
const SOS = SwitchOnSafety;

using MathOptInterface
const MOI = MathOptInterface

function symsolve(M, set_variable; volume_heuristic=set -> L1_heuristic(set, sym_rect(M)))
    hs = system(M, T, N=1, sym=true, vmax = 1.0, v=(1.0,))
    # Shift interval 5 <= v <= 15.6 -> -5.3 <= v <= 5.3
    sets = invariant_sets(hs, sdp_factory, [set_variable], volume_heuristic=volume_heuristic)
    [SOS.SetProg.Sets.Translation(set, _vec(M, 0, 2.0, 0)) for set in sets]
end

function fullsolve(T, M, setvar::Function;
                   volume_heuristic=set -> L1_heuristic(set, sym_rect(M)),
                   oneshot = false,
                   onlyone = false)
    hs = system(M, T)
    function hv(v)
        h = zeros(statedim(hs, 1))
        for i in 1:M
            h[2i] = (vmin .+ v) / 2
        end
        h[2M+1] = (vmin .+ v) / 2
        h
    end
    habc = SOS.SetProg.InteriorPoint.(hv.(v))
    ha = habc[1]
    hi = [ha, ha, ha, ha, ha, ha, ha, ha, ha, ha]
    set_variables = map(setvar, hi)
    if oneshot
        λ = Dict(t => 1.0 for t in transitions(hs))
    else
        λ = Dict(t => 1.0 for t in transitions(hs) if source(hs, t) == target(hs, t))
    end
    if oneshot
        return invariant_sets(hs, sdp_factory, set_variables, λ = λ, volume_heuristic=volume_heuristic);
    else
        sets = Vector{SwitchOnSafety.SetProg.Sets.AbstractSet{Float64}}(undef, N)
        for i in (onlyone ? (1:1) : (1:N))
            invariant_sets!(sets, i:i, hs, sdp_factory, set_variables, λ = λ, enabled=1:i,
                            volume_heuristic=volume_heuristic)
        end
        return sets
    end
end

using HybridSystems
include(joinpath(dirname(dirname(pathof(HybridSystems))), "examples", "cruise_control.jl"));

if true
    D = 1.0
    U = 1.0
    v_shift = 2.0
    vmin = -1.0
    vmax = 2.0
    v = (1.0,)
    m = 1.0
    m0 = 1.0
    h = 0.5
    kd = 1/2
    ks = 1/2
    Δv = 5.0
else
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
    vmax = 35.0
    m = 1000
    m0 = 500
    h = 0.8
    kd = 4600
    ks = 4500
    Δv = (v[1] - vmin) / 2
end

T = 2
N = 1 + (T+1) * length(v)
function system(M, T; v=v, N = 1 + (T+1) * length(v), sym=false, vmax = vmax)
    H = h * T
    #return cruise_control_example(N, M, vmin = vmin, vmax=vmax, v=v, U=U, H=H, D=D, T=T, sym=sym, m0 = m0)
    return cruise_control_example(N, M, vmin = vmin, vmax=vmax, v=v,
        U=U, H=h*T, D=D, T=T, sym=sym, m0 = m0, ks=ks, m=m, kd=kd)
end
_vec(M, d, v, u) = [repeat([d, v], M); v; u]
# |string_elongation| < D
# |acceleration| < U
sym_rect(M) = _vec(M, D, 1.0, U) # _vec(M, D, Δv, U)

using SetProg
using SwitchOnSafety
const SOS = SwitchOnSafety;

using MathOptInterface
const MOI = MathOptInterface

function symsolve(M, set_variable; volume_heuristic=set -> L1_heuristic(set, sym_rect(M)))
    hs = system(M, T, N=1, sym=true)
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
    @show nstates(hs)
    hi = fill(ha, nstates(hs))
    set_variables = map(setvar, hi)
    if oneshot
        λ = Dict(t => 1.0 for t in transitions(hs))
    else
        λ = Dict(t => 1.0 for t in transitions(hs) if source(hs, t) == target(hs, t))
    end
    if oneshot
        return invariant_sets(hs, sdp_factory, set_variables, λ = λ, volume_heuristic=volume_heuristic);
    else
        sets = Vector{SwitchOnSafety.SetProg.Sets.AbstractSet{Float64}}(undef, nstates(hs))
        for i in (onlyone ? (1:1) : (1:nstates(hs)))
            println("Computing set $i")
            invariant_sets!(sets, i:i, hs, sdp_factory, set_variables, λ = λ, enabled=1:i,
                            volume_heuristic=volume_heuristic)
        end
        return sets
    end
end

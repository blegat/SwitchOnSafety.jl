import GLPK
solver = GLPK.Optimizer;
import QHull
Qlib = QHull.Library(solver);
import SwitchOnSafety
include(joinpath(dirname(dirname(pathof(SwitchOnSafety))), "examples", "cruise_control.jl"));

function fixed_point_iteration(set::Polyhedron, sys::LinearAlgebraicDiscreteSystem)
    new_set = set ∩ (sys.A \ (sys.E * set))
    removehredundancy!(new_set)
    return new_set::typeof(set)
end
iterate!(sets, sys) = push!(sets, fixed_point_iteration(sets[end], sys))
function iterate!(sets, sys, n)
    if n > 0
        iterate!(sets, sys)
        iterate!(sets, sys, n - 1)
    end
end

function _proj(P, I)
    Q = project(P, I)
    removehredundancy!(Q)
    return Q
end
_proj(Ps::Vector, I) = map(P -> _proj(P, I), Ps)

hsys = cruise_control_example(1, 1, vmin=vmin, vmax=vmax, v=v,
                              U=U, H=1/4, D=D, T=T, sym=true, m0 = m0, ks=1, m=2m0, kd=1/4, lib=Qlib)
sys_control = hsys.resetmaps[1]
sys_alg = SwitchOnSafety.algebraiclift(sys_control);

QX = [hsys.modes[1].X]
iterate!(QX, sys_alg, 10)
QX123 = _proj(QX, 1:3)
nhalfspaces.(QX123)

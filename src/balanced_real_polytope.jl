using JuMP
using ParameterJuMP
using Polyhedra
mutable struct BalancedRealPolytope{T, VT <: AbstractVector{T}, D<:Polyhedra.FullDim} <: Polyhedra.VPolytope{T}
    d::D
    points::Vector{VT}
    optimizer_constructor
    model::Union{Nothing, MOI.ModelLike}
    sum_con::Union{Nothing, MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, MOI.LessThan{T}}}
    z::Union{Nothing, Vector{MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, MOI.EqualTo{T}}}}
    t_0::Union{Nothing, MOI.VariableIndex}
    function BalancedRealPolytope{T, VT, D}(
        d::Polyhedra.FullDim, points::Polyhedra.PointIt,
        optimizer_constructor=nothing) where {T, VT, D}
        new{T, VT, D}(Polyhedra.FullDim_convert(D, d), Polyhedra.lazy_collect(points), optimizer_constructor, nothing, nothing, nothing, nothing)
    end
end
function BalancedRealPolytope(d::Polyhedra.FullDim, points::Polyhedra.PointIt, args...)
    return BalancedRealPolytope{Polyhedra.coefficient_type(eltype(points)), eltype(points),
                                typeof(d)}(d, points, args...)
end

Polyhedra.FullDim(v::BalancedRealPolytope) = v.d
Polyhedra.vvectortype(::Type{<:BalancedRealPolytope{T, AT}}) where {T, AT} = AT
function Polyhedra.similar_type(PT::Type{<:BalancedRealPolytope}, d::Polyhedra.FullDim, ::Type{T}) where {T}
    BalancedRealPolytope{T, Polyhedra.similar_type(Polyhedra.vvectortype(PT), d, T), typeof(d)}
end

Polyhedra.vreptype(::Type{BalancedRealPolytope{T, AT, D}}) where {T, AT, D} = Polyhedra.Hull{T, AT, D}
Polyhedra.fulltype(::Type{BalancedRealPolytope{T, AT, D}}) where {T, AT, D} = Polyhedra.Hull{T, AT, D}

Base.length(idxs::Polyhedra.PointIndices{T, <:BalancedRealPolytope{T}}) where {T} = 2length(idxs.rep.points)
Base.isempty(idxs::Polyhedra.PointIndices{T, <:BalancedRealPolytope{T}}) where {T} = isempty(idxs.rep.points)
function Polyhedra.startindex(idxs::Polyhedra.PointIndices{T, <:BalancedRealPolytope{T}}) where {T}
    if isempty(idxs.rep.points)
        return nothing
    else
        return eltype(idxs)(1)
    end
end
function Base.get(rep::BalancedRealPolytope{T}, idx::Polyhedra.PointIndex{T}) where {T}
    return sign(idx.value) * rep.points[abs(idx.value)]
end
function Polyhedra.nextindex(rep::BalancedRealPolytope{T}, idx::Polyhedra.PointIndex{T}) where {T}
    if idx.value < 0
        if -idx.value >= length(rep.points)
            return nothing
        else
            return typeof(idx)(idx.value - 1)
        end
    else
        if idx.value >= length(rep.points)
            return typeof(idx)(-1)
        else
            return typeof(idx)(idx.value + 1)
        end
    end
end

function _abs_constraint(model::MOI.ModelLike, t::MOI.VariableIndex, q::MOI.VariableIndex, T::Type)
    MOI.add_constraint(model, MOI.Utilities.operate(-, T, q, t), MOI.GreaterThan(zero(T)))
    MOI.add_constraint(model, MOI.Utilities.operate(+, T, q, t), MOI.GreaterThan(zero(T)))
end

function _build_model(brp::BalancedRealPolytope{T}) where T
    @assert brp.model === nothing
    n = length(brp.points)
    brp.model = MOI.instantiate(brp.optimizer_constructor, with_bridge_type = T)
    t = MOI.add_variables(brp.model, n)
    q = MOI.add_variables(brp.model, n)
    # |t| ≤ q
    for i in 1:n
        _abs_constraint(brp.model, t[i], q[i], T)
    end
    # sum |t| ≤ sum q ≤ 1
    brp.sum_con = MOI.add_constraint(brp.model, _sum(q, T), MOI.LessThan(one(T)))
    return _convex_combination(t, brp.points, brp.d)
end

function _parametrized_model(brp::BalancedRealPolytope, v)
    hull = _build_model(brp)
    brp.z = [MOI.add_constraint(brp.model, hull[i], MOI.EqualTo(v[i])) for i in 1:brp.d]
end
function _ratio_model(brp::BalancedRealPolytope{T}, v) where T
    hull = _build_model(brp)
    brp.t_0 = MOI.add_variable(brp.model)
    for i in 1:brp.d
        push!(hull[i].terms, MOI.ScalarAffineTerm(-v[i], brp.t_0))
    end
    brp.z = [MOI.add_constraint(brp.model, hull[i], MOI.EqualTo(zero(T))) for i in 1:brp.d]
end
function _update_model(p::BalancedRealPolytope{T}, v) where T
    t = MOI.add_variable(p.model)
    q = MOI.add_variable(p.model)
    _abs_constraint(p.model, t, q, T)
    MOI.modify(p.model, p.sum_con, MOI.ScalarCoefficientChange(q, one(T)))
    for i in 1:p.d
        MOI.modify(p.model, p.z[i], MOI.ScalarCoefficientChange(t, v[i]))
    end
end

function _fix(brp::BalancedRealPolytope, v)
    for i in 1:brp.d
        if brp.t_0 === nothing
            MOI.set(brp.model, MOI.ConstraintSet(), brp.z[i], MOI.EqualTo(v[i]))
        else
            MOI.modify(brp.model, brp.z[i], MOI.ScalarCoefficientChange(brp.t_0, -v[i]))
        end
    end
end

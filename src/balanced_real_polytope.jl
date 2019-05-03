using JuMP
using ParameterJuMP
using Polyhedra
mutable struct BalancedRealPolytope{T, VT <: AbstractVector{T}, D<:Polyhedra.FullDim} <: Polyhedra.VPolytope{T}
    d::D
    points::Vector{VT}
    factory::Union{Nothing, JuMP.OptimizerFactory}
    model::Union{Nothing, JuMP.Model}
    z::Union{Nothing, Vector{ParameterJuMP.ParameterRef}}
    t_0::Union{Nothing, JuMP.VariableRef}
    function BalancedRealPolytope{T, VT, D}(
        d::Polyhedra.FullDim, points::Polyhedra.PointIt,
        factory::Union{Nothing, JuMP.OptimizerFactory}=nothing) where {T, VT, D}
        new{T, VT, D}(Polyhedra.FullDim_convert(D, d), Polyhedra.lazy_collect(points), factory, nothing, nothing, nothing)
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

function _build_model(brp::BalancedRealPolytope)
    n = length(brp.points)
    brp.model = ParameterJuMP.ModelWithParams(brp.factory)
    t = @variable(brp.model, [1:n])
    q = @variable(brp.model, [1:n])
    # |t| ≤ q
    @constraint(brp.model,  t .≤ q)
    @constraint(brp.model, -t .≤ q)
    # sum |t| ≤ sum q ≤ 1
    @constraint(brp.model, sum(q) ≤ 1)
    return sum(t[i] .* brp.points[i] for i in 1:n)
end

function _parametrized_model(brp::BalancedRealPolytope, v)
    hull = _build_model(brp)
    brp.z = ParameterJuMP.add_parameters(brp.model, collect(v))
    @constraint(brp.model, brp.z .== hull)
end

function _ratio_model(brp::BalancedRealPolytope, v)
    hull = _build_model(brp)
    brp.t_0 = @variable(brp.model)
    @constraint(brp.model, brp.t_0 .* v .== hull)
end

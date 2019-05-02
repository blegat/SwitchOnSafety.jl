using JuMP
using ParameterJuMP
using Polyhedra

mutable struct Conitope{T, VT <: AbstractVector{T}, D<:Polyhedra.FullDim}
    d::D
    points::Vector{VT}
    factory::Union{Nothing, JuMP.OptimizerFactory}
    model::Union{Nothing, JuMP.Model}
    z::Union{Nothing, Vector{ParameterJuMP.ParameterRef}}
    t_0::Union{Nothing, JuMP.VariableRef}
    function Conitope{T, VT, D}(
        d::Polyhedra.FullDim, points::Polyhedra.PointIt,
        factory::Union{Nothing, JuMP.OptimizerFactory}=nothing) where {T, VT, D}
        new{T, VT, D}(Polyhedra.FullDim_convert(D, d),
                      Polyhedra.lazy_collect(points), factory, nothing, nothing, nothing)
    end
end
function Conitope(d::Polyhedra.FullDim, points::Polyhedra.PointIt, args...)
    return Conitope{Polyhedra.coefficient_type(eltype(points)), eltype(points),
                    typeof(d)}(d, points, args...)
end

Polyhedra.FullDim(v::Conitope) = v.d

function _build_model(brp::Conitope)
    n = length(brp.points)
    brp.model = ParameterJuMP.ModelWithParams(brp.factory)
    @variable(brp.model, t[1:n] ≥ 0)
    @constraint(brp.model, sum(t) ≤ 1)
    return sum(t[i] .* brp.points[i] for i in 1:n)
end
function _parametrized_model(brp::Conitope, v)
    hull = _build_model(brp)
    brp.z = ParameterJuMP.add_parameters(brp.model, collect(v))
    @constraint(brp.model, brp.z .≤ hull)
end
function _ratio_model(brp::Conitope, v)
    hull = _build_model(brp)
    brp.t_0 = @variable(brp.model)
    @constraint(brp.model, brp.t_0 .* v .≤ hull)
end

function _ei(i::Integer, d::Integer, T::Type)
    e = zeros(T, d)
    e[i] = one(T)
    return e
end
function Polyhedra.polyhedron(
    c::Conitope{T},
    lib::Polyhedra.Library=Polyhedra.default_library(Polyhedra.FullDim(c),
                                                     T)) where T
    n = fulldim(c)
    p = polyhedron(vrep([Ray(-_ei(i, n, T)) for i in 1:n]) + vrep(c.points),
                   lib)
    return p ∩ (Matrix(-LinearAlgebra.I, n, n) * p)
end

using JuMP
using ParameterJuMP
using Polyhedra
mutable struct BalancedComplexPolytope{T, VT <: AbstractVector{T}, D<:Polyhedra.FullDim}
    d::D
    points::Vector{VT}
    factory::Union{Nothing, JuMP.OptimizerFactory}
    model::Union{Nothing, JuMP.Model}
    real_z::Union{Nothing, Vector{ParameterJuMP.ParameterRef}}
    imag_z::Union{Nothing, Vector{ParameterJuMP.ParameterRef}}
    t_0::Union{Nothing, JuMP.VariableRef}
    function BalancedComplexPolytope{T, VT, D}(
        d::Polyhedra.FullDim, points::Polyhedra.PointIt,
        factory::Union{Nothing, JuMP.OptimizerFactory}=nothing) where {T, VT, D}
        new{T, VT, D}(Polyhedra.FullDim_convert(D, d),
                      Polyhedra.lazy_collect(points), factory, nothing, nothing, nothing, nothing)
    end
end
function BalancedComplexPolytope(d::Polyhedra.FullDim, points::Polyhedra.PointIt, args...)
    return BalancedComplexPolytope{Polyhedra.coefficient_type(eltype(points)), eltype(points),
                                   typeof(d)}(d, points, args...)
end

Polyhedra.FullDim(v::BalancedComplexPolytope) = v.d

function _reset_model(p::BalancedComplexPolytope)
    p.model = nothing
    p.real_z = nothing
    p.imag_z = nothing
    p.t_0 = nothing
end

function _build_model(brp::BalancedComplexPolytope)
    n = length(brp.points)
    brp.model = ParameterJuMP.ModelWithParams(brp.factory)
    # Real part of `t`
    α = @variable(brp.model, [1:n])
    # Imaginary part of `t`
    β = @variable(brp.model, [1:n])
    q = @variable(brp.model, [1:n])
    # |t| ≤ q
    for i in 1:n
        @constraint(brp.model,  [q[i], α[i], β[i]] in SecondOrderCone())
    end
    # sum |t| ≤ sum q ≤ 1
    @constraint(brp.model, sum(q) ≤ 1)
    return sum(α[i] .* real.(brp.points[i]) for i in 1:n) - sum(β[i] .* imag.(brp.points[i]) for i in 1:n),
           sum(α[i] .* imag.(brp.points[i]) for i in 1:n) + sum(β[i] .* real.(brp.points[i]) for i in 1:n)
end

function _fix(p::BalancedComplexPolytope, v)
    JuMP.fix.(brp.real_z, real.(v))
    JuMP.fix.(brp.imag_z, imag.(v))
end

function _parametrized_model(brp::BalancedComplexPolytope, v)
    real_hull, imag_hull = _build_model(brp)
    brp.real_z = ParameterJuMP.add_parameters(brp.model, collect(real.(v)))
    brp.imag_z = ParameterJuMP.add_parameters(brp.model, collect(imag.(v)))
    @constraint(brp.model, brp.real_z .== real_hull)
    @constraint(brp.model, brp.imag_z .== imag_hull)
end

function _ratio_model(brp::BalancedComplexPolytope, v)
    real_hull, imag_hull = _build_model(brp)
    brp.t_0 = @variable(brp.model)
    @constraint(brp.model, brp.t_0 .* real.(v) .== real_hull)
    @constraint(brp.model, brp.t_0 .* imag.(v) .== imag_hull)
end

using JuMP
using ParameterJuMP
using Polyhedra

mutable struct Conitope{T, VT <: AbstractVector{T}, D<:Polyhedra.FullDim}
    d::D
    points::Vector{VT}
    optimizer_constructor
    model::Union{Nothing, MOI.ModelLike}
    sum_con::Union{Nothing, MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, MOI.LessThan{T}}}
    z::Union{Nothing, Vector{MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, MOI.GreaterThan{T}}}}
    t_0::Union{Nothing, MOI.VariableIndex}
    function Conitope{T, VT, D}(
        d::Polyhedra.FullDim, points::Polyhedra.PointIt,
        optimizer_constructor=nothing) where {T, VT, D}
        new{T, VT, D}(Polyhedra.FullDim_convert(D, d),
                      Polyhedra.lazy_collect(points), optimizer_constructor, nothing, nothing, nothing)
    end
end
function Conitope(d::Polyhedra.FullDim, points::Polyhedra.PointIt, args...)
    return Conitope{Polyhedra.coefficient_type(eltype(points)), eltype(points),
                    typeof(d)}(d, points, args...)
end

Polyhedra.FullDim(v::Conitope) = v.d

function _sum(x::Vector{MOI.VariableIndex}, T::Type)
    MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(one(T), xi) for xi in x],
                             zero(T))
end
function _convex_combination(t::Vector{MOI.VariableIndex}, points::Vector{Vector{T}}, d::Int) where T
    return [MOI.ScalarAffineFunction([
        MOI.ScalarAffineTerm(points[i][j], t[i]) for i in eachindex(t)
    ], zero(T)) for j in 1:d]
end
function _build_model(brp::Conitope{T}) where T
    @assert brp.model === nothing
    n = length(brp.points)
    brp.model = _instantiate(brp.optimizer_constructor, T)
    t, ct = MOI.add_constrained_variables(brp.model, MOI.Nonnegatives(n))
    # sum t_i ≤ 1
    brp.sum_con = MOI.add_constraint(brp.model, _sum(t, T), MOI.LessThan(one(T)))
    return _convex_combination(t, brp.points, brp.d)
end
function _parametrized_model(brp::Conitope, v)
    hull = _build_model(brp)
    brp.z = [MOI.add_constraint(brp.model, hull[i], MOI.GreaterThan(v[i])) for i in 1:brp.d]
end
function _ratio_model(brp::Conitope{T}, v) where T
    hull = _build_model(brp)
    brp.t_0 = MOI.add_variable(brp.model)
    for i in 1:brp.d
        push!(hull[i].terms, MOI.ScalarAffineTerm(-v[i], brp.t_0))
    end
    brp.z = [MOI.add_constraint(brp.model, hull[i], MOI.GreaterThan(zero(T))) for i in 1:brp.d]
end
function _update_model(p::Conitope{T}, v) where T
    ts, ct = MOI.add_constrained_variables(p.model, MOI.Nonnegatives(1))
    t = ts[1]
    MOI.modify(p.model, p.sum_con, MOI.ScalarCoefficientChange(t, one(T)))
    for i in 1:p.d
        MOI.modify(p.model, p.z[i], MOI.ScalarCoefficientChange(t, v[i]))
    end
end

function _fix(brp::Conitope, v)
    for i in 1:brp.d
        if brp.t_0 === nothing
            MOI.set(brp.model, MOI.ConstraintSet(), brp.z[i], MOI.GreaterThan(v[i]))
        else
            MOI.modify(brp.model, brp.z[i], MOI.ScalarCoefficientChange(brp.t_0, -v[i]))
        end
    end
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

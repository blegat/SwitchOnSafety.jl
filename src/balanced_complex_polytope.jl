using JuMP
using ParameterJuMP
using Polyhedra
mutable struct BalancedComplexPolytope{T, CT, VT <: AbstractVector{CT}, D<:Polyhedra.FullDim}
    d::D
    points::Vector{VT}
    optimizer_constructor
    model::Union{Nothing, MOI.ModelLike}
    sum_con::Union{Nothing, MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, MOI.LessThan{T}}}
    real_z::Union{Nothing, Vector{MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, MOI.EqualTo{T}}}}
    imag_z::Union{Nothing, Vector{MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, MOI.EqualTo{T}}}}
    t_0::Union{Nothing, MOI.VariableIndex}
    function BalancedComplexPolytope{T, CT, VT, D}(
        d::Polyhedra.FullDim, points::Polyhedra.PointIt,
        optimizer_constructor=nothing) where {T, CT, VT, D}
        new{T, CT, VT, D}(Polyhedra.FullDim_convert(D, d),
                          Polyhedra.lazy_collect(points), optimizer_constructor, nothing, nothing, nothing, nothing)
    end
end
function BalancedComplexPolytope(d::Polyhedra.FullDim, points::Polyhedra.PointIt, args...)
    VT = eltype(points)
    CT = Polyhedra.coefficient_type(VT)
    return BalancedComplexPolytope{real(CT), CT, VT, typeof(d)}(d, points, args...)
end

Polyhedra.FullDim(v::BalancedComplexPolytope) = v.d

function _build_model(brp::BalancedComplexPolytope{T}) where T
    @assert brp.model === nothing
    n = length(brp.points)
    brp.model = MOI.instantiate(brp.optimizer_constructor, with_bridge_type = T)
    # Real part of `t`
    α = Vector{MOI.VariableIndex}(undef, n)
    # Imaginary part of `t`
    β = Vector{MOI.VariableIndex}(undef, n)
    q = Vector{MOI.VariableIndex}(undef, n)
    # |t| ≤ q
    for i in 1:n
        vars, con = MOI.add_constrained_variables(brp.model, MOI.SecondOrderCone(3))
        q[i], α[i], β[i] = vars
    end
    # sum |t| ≤ sum q ≤ 1
    brp.sum_con = MOI.add_constraint(brp.model, _sum(q, T), MOI.LessThan(one(T)))
    rp = [real.(p) for p in brp.points]
    ip = [imag.(p) for p in brp.points]
    real_hull = MOI.Utilities.operate!.(
        -, T,
        _convex_combination(α, rp, brp.d),
        _convex_combination(β, ip, brp.d))
    imag_hull = MOI.Utilities.operate!.(
        +, T,
        _convex_combination(α, ip, brp.d),
        _convex_combination(β, rp, brp.d))
    return real_hull, imag_hull
end

function _parametrized_model(brp::BalancedComplexPolytope, v)
    real_hull, imag_hull = _build_model(brp)
    brp.real_z = [MOI.add_constraint(brp.model, real_hull[i], MOI.EqualTo(real(v[i]))) for i in 1:brp.d]
    brp.imag_z = [MOI.add_constraint(brp.model, imag_hull[i], MOI.EqualTo(imag(v[i]))) for i in 1:brp.d]
end

function _ratio_model(brp::BalancedComplexPolytope{T}, v) where T
    real_hull, imag_hull = _build_model(brp)
    brp.t_0 = MOI.add_variable(brp.model)
    for i in 1:brp.d
        push!(real_hull[i].terms, MOI.ScalarAffineTerm(-real(v[i]), brp.t_0))
        push!(imag_hull[i].terms, MOI.ScalarAffineTerm(-imag(v[i]), brp.t_0))
    end
    brp.real_z = [MOI.add_constraint(brp.model, real_hull[i], MOI.EqualTo(zero(T))) for i in 1:brp.d]
    brp.imag_z = [MOI.add_constraint(brp.model, imag_hull[i], MOI.EqualTo(zero(T))) for i in 1:brp.d]
end

function _update_model(p::BalancedComplexPolytope{T}, v) where T
    vars, con = MOI.add_constrained_variables(p.model, MOI.SecondOrderCone(3))
    q, α, β = vars
    MOI.modify(p.model, p.sum_con, MOI.ScalarCoefficientChange(q, one(T)))
    for i in 1:p.d
        MOI.modify(p.model, p.real_z[i], MOI.ScalarCoefficientChange(α, real(v[i])))
        MOI.modify(p.model, p.imag_z[i], MOI.ScalarCoefficientChange(β, imag(v[i])))
    end
end

function _fix(bcp::BalancedComplexPolytope, v)
    for i in 1:bcp.d
        if bcp.t_0 === nothing
            MOI.set(bcp.model, MOI.ConstraintSet(), bcp.real_z[i], MOI.EqualTo(real(v[i])))
            MOI.set(bcp.model, MOI.ConstraintSet(), bcp.imag_z[i], MOI.EqualTo(imag(v[i])))
        else
            MOI.modify(bcp.model, bcp.real_z[i], MOI.ScalarCoefficientChange(bcp.t_0, -real(v[i])))
            MOI.modify(bcp.model, bcp.imag_z[i], MOI.ScalarCoefficientChange(bcp.t_0, -imag(v[i])))
        end
    end
end

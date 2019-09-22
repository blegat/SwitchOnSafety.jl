export invariant_set, invariant_sets, invariant_sets!, algebraiclift

using FillArrays

#function r(A::Matrix{T}, c::Vector{T}=zeros(T, size(A, 1))) where T
#    [one(T) zeros(T, 1, size(A, 2))
#     c A]
#end

# x+ = Ax + Bu, x ∈ X, u ∈ U
# y+ = [A B; 0 0]y [0; I]w, y ∈ X × U
# [I 0]y+ = [A B]y, y ∈ X × U
function algebraiclift(s::ConstrainedLinearControlDiscreteSystem)
    A = [s.A s.B]
    E = [Matrix(one(eltype(s.A)) * I, size(s.A)...) zeros(eltype(s.A), size(s.B)...)]
    return ConstrainedLinearAlgebraicDiscreteSystem(A, E, stateset(s) * inputset(s))
end
function algebraiclift(s::ConstrainedLinearControlDiscreteSystem{T, MTA, MTB, ST, FullSpace}) where {T, MTA <: AbstractMatrix{T}, MTB <: AbstractMatrix{T}, ST}
    as = algebraiclift(LinearControlDiscreteSystem(s.A, s.B))
    return ConstrainedLinearAlgebraicDiscreteSystem(as.A, as.E, stateset(s))
end
function algebraiclift(s::LinearControlDiscreteSystem)
    n = statedim(s)
    z = findall(i -> iszero(sum(abs.(s.B[i,:]))), 1:n)
    # TODO ty - 1//2y^3 + 3//1xy + 2//1yhe affine space may not be parallel to classical axis
    LinearAlgebraicDiscreteSystem(s.A[z, :], Matrix(one(eltype(s.A)) * I, n, n)[z, :])
end
algebraiclift(s::ConstrainedDiscreteIdentitySystem) = s
algebraiclift(S::AbstractVector) = algebraiclift.(S)
algebraiclift(S::Fill) = Fill(algebraiclift(first(S)), length(S))
function algebraiclift(h::HybridSystem)
    HybridSystem(h.automaton, algebraiclift(h.modes), algebraiclift(h.resetmaps), h.switchings)
end

const DTAHAS = HybridSystem{<:AbstractAutomaton, <:ConstrainedContinuousIdentitySystem, <:LinearAlgebraicDiscreteSystem}
const DTAHCS = HybridSystem{<:AbstractAutomaton, <:ConstrainedContinuousIdentitySystem, <:LinearControlDiscreteSystem}

function constrain_invariance(model, ::ContinuousIdentitySystem, set) end
function constrain_invariance(model, s::ConstrainedContinuousIdentitySystem, set)
    return @constraint(model, set ⊆ stateset(s))
end

function constrain_invariance(model, ::IdentityMap, source_set, target_set, λ) end
function constrain_invariance(model, s::LinearMap, source_set, target_set, λ)
    return @constraint(model, s.A * source_set ⊆ target_set)
end
function constrain_invariance(model, s::LinearAlgebraicDiscreteSystem, source_set, target_set, λ)
    return @constraint(model, s.A * source_set ⊆ s.E * target_set,
                       S_procedure_scaling = λ)
end

function isinfeasible(status::Tuple{MOI.TerminationStatusCode, MOI.ResultStatusCode, MOI.ResultStatusCode})
    status[3] == MOI.INFEASIBILITY_CERTIFICATE
end
function isfeasible(status::Tuple{MOI.TerminationStatusCode, MOI.ResultStatusCode, MOI.ResultStatusCode})
    status[2] == MOI.FEASIBLE_POINT
end
function isdecided(status::Tuple{MOI.TerminationStatusCode, MOI.ResultStatusCode, MOI.ResultStatusCode})
    return isinfeasible(status) || isfeasible(status)
end

"""
    invariant_sets!(sets, modes_to_compute, system::AbstractHybridSystem,
                    args...; kwargs...)

Similar to `invariant_sets(system, args...; kwargs...)` but stores the result
in `sets` and only compute the modes in `modes_to_compute`. The other sets in
`enabled` that are not in `modes_to_compute` are assumed to have the value given
in `sets`.
"""
function invariant_sets! end

"""
    invariant_sets(system::AbstractHybridSystem, factory::JuMP.OptimizerFactory,
                   set_variables::AbstractVector{<:SetProg.AbstractVariable};
                   volume_heuristic = nth_root,
                   infeasibility_certificates = nothing,
                   verbose=1,
                   λ=Dict{transitiontype(system), Float64}(),
                   enabled = states(system))

Compute maximal invariant sets of the family `set_variables` for the modes of
`system` using the solver provided by `factory`. The volume of the sets is
estimated using `volume_heuristic`. If the program is infeasible, the
certificates for each transition are stored in `infeasibility_certificates`.
For the containment of non-homogeneous, the S-procedure might be a Bilinear
Matrix Inequality (BMI) which is NP-hard to solve. To avoid that, provides the
value of `λ` to use in the dictionary `λ`. To ignore some state and the
transitions involving these states in the computation, give an `enabled` vector
without them.
"""
function invariant_sets end


function invariant_sets!(
    sets, modes_to_compute, s::AbstractHybridSystem, factory::JuMP.OptimizerFactory,
    set_variables::AbstractVector{<:SetProg.AbstractVariable} = map(
        cv->Ellipsoid(point=SetProg.InteriorPoint(cv[1])),
        chebyshevcenter.(stateset.(s.modes)));
    λ=Dict{transitiontype(s), Float64}(),
    enabled = states(s),
    volume_heuristic = nth_root,
    infeasibility_certificates = nothing,
    verbose=1
)
    model = SOSModel(factory)
    #set_vrefs = @variable(model, [q in modes_to_compute], Ellipsoid(point=h[q]))
    # FIXME calls JuMP.variable_type(model, Ellipsoid(point=h[q])) then
    # complains that q is not defined
    set_vrefs = Dict(state => @variable(model, variable_type=set_variables[state]) for state in modes_to_compute)

    if volume_heuristic !== nothing
        @objective(model, Max,
                   sum(set -> volume_heuristic(volume(set)), values(set_vrefs)))
    end

    λouts = Dict{transitiontype(s), JuMP.AffExpr}()
    transition_constraints = HybridSystems.transition_property(
        s, ConstraintRef{Model, SetProg.ConstraintIndex, SetProg.SetShape})

    for q in modes_to_compute
        source_set = set_vrefs[q]
        # Invariance constraint for mode
        constrain_invariance(model, HybridSystems.mode(s, q), source_set)
        # Invariance constraint for transitions
        for t in out_transitions(s, q)
            λin = get(λ, t, nothing)
            target_state = target(s, t)
            if target_state in enabled
                target_set = target_state in modes_to_compute ? set_vrefs[target_state] : sets[target_state]
                if λin !== nothing
                    λouts[t] = λin
                end
                con_ref = constrain_invariance(
                    model, resetmap(s, t), source_set, target_set, λin)
                if con_ref !== nothing
                    transition_constraints[t] = con_ref
                end
            end
        end
    end

    if verbose >= 3
        print(model)
    end

    JuMP.optimize!(model)

    if verbose >= 4
        print(model)
    end

    if verbose >= 1
        try
            @show MOI.get(model, MOI.SolveTime())
        catch err
            if !(err isa ArgumentError)
                rethrow(err)
            end
        end
        @show JuMP.termination_status(model)
        @show JuMP.primal_status(model)
        @show JuMP.dual_status(model)
        @show JuMP.objective_value(model)
    end

    if verbose >= 2
        for (t, λout) in λouts
            println("λ for $t is $(JuMP.value(λout))")
        end
    end

    status = (JuMP.termination_status(model),
              JuMP.primal_status(model),
              JuMP.dual_status(model))

    if isfeasible(status)
        for q in modes_to_compute
            sets[q] = JuMP.value(set_vrefs[q])
        end
    elseif isinfeasible(status) && infeasibility_certificates !== nothing
        for q in modes_to_compute
            for t in out_transitions(s, q)
                if target(s, t) in enabled
                    infeasibility_certificates[t] = SetProg.SumOfSquares.moment_matrix(transition_constraints[t])
                end
            end
        end
    end

    return status
end
function invariant_sets!(sets, modes_to_compute, s::DTAHCS, args...; kws...)
    return invariant_sets!(sets, modes_to_compute, algebraiclift(s), args...;
                           kws...)
end

function invariant_sets(s::AbstractHybridSystem, args...; kws...)
    sets = HybridSystems.state_property(s, SetProg.Sets.AbstractSet{Float64})
    invariant_sets!(sets, states(s), s, args...; kws...)
    return sets
end
function invariant_sets(s::DTAHCS, args...; kws...)
    return invariant_sets(algebraiclift(s), args...; kws...)
end

_mode(s::AbstractContinuousSystem) = s
_mode(s::AbstractDiscreteSystem) = ConstrainedContinuousIdentitySystem(statedim(s), stateset(s))

_resetmap(s::AbstractContinuousSystem) = IdentityMap(statedim(s))
_resetmap(s::ConstrainedDiscreteIdentitySystem) = IdentityMap(statedim(s))
_resetmap(s::ConstrainedLinearDiscreteSystem) = LinearMap(s.A)
# TODO need to add LinearAlgebraicMap to MathematicalSystems
_resetmap(s::ConstrainedLinearAlgebraicDiscreteSystem) = LinearAlgebraicDiscreteSystem(s.A, s.E)

function invariant_set(
    system::Union{
        ConstrainedContinuousIdentitySystem,
        ConstrainedDiscreteIdentitySystem,
        ConstrainedLinearDiscreteSystem,
        ConstrainedLinearAlgebraicDiscreteSystem},
    factory::JuMP.OptimizerFactory,
    set_variable::SetProg.AbstractVariable;
    λ=nothing,
    kws...)
    hs = HybridSystem(
        HybridSystems.OneStateAutomaton(1),
        [_mode(system)], [_resetmap(system)], [AutonomousSwitching()])
    λs = Dict{transitiontype(hs), Float64}()
    if λ != nothing
        λs[HybridSystems.OneStateTransition(1)] = λ
    end
    sets = invariant_sets(hs, factory, [set_variable]; λ = λs, kws...)
    return sets[1]
end
function invariant_set(s::ConstrainedLinearControlDiscreteSystem{T, MTA, MTB, ST, FullSpace}, args...; kws...) where {T, MTA <: AbstractMatrix{T}, MTB <: AbstractMatrix{T}, ST}
    return invariant_set(algebraiclift(s), args...; kws...)
end
function invariant_set(s::ConstrainedLinearControlDiscreteSystem, args...; kws...)
    set = invariant_set(algebraiclift(s), args...; kws...)
    return project(set, 1:statedim(s))
end

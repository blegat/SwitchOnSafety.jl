export invariant_sets, invariant_sets!, algebraiclift

using FillArrays

#function r(A::Matrix{T}, c::Vector{T}=zeros(T, size(A, 1))) where T
#    [one(T) zeros(T, 1, size(A, 2))
#     c A]
#end

function algebraiclift(s::LinearControlDiscreteSystem)
    n = statedim(s)
    z = findall(i -> iszero(sum(abs.(s.B[i,:]))), 1:n)
    # TODO ty - 1//2y^3 + 3//1xy + 2//1yhe affine space may not be parallel to classical axis
    LinearAlgebraicDiscreteSystem(s.A[z, :], Matrix(1.0I, n, n)[z, :])
end
algebraiclift(s::ConstrainedDiscreteIdentitySystem) = s
algebraiclift(S::AbstractVector) = algebraiclift.(S)
algebraiclift(S::Fill) = Fill(algebraiclift(first(S)), length(S))
function algebraiclift(h::HybridSystem)
    HybridSystem(h.automaton, algebraiclift(h.modes), algebraiclift(h.resetmaps), h.switchings)
end

const DTAHAS = HybridSystem{<:AbstractAutomaton, <:ConstrainedDiscreteIdentitySystem, <:LinearAlgebraicDiscreteSystem}
const DTAHCS = HybridSystem{<:AbstractAutomaton, <:ConstrainedDiscreteIdentitySystem, <:LinearControlDiscreteSystem}
function _vars(s::DTAHAS)
    @polyvar x[1:statedim(s, 1)] z
    [z; x]
end

function invariant_sets!(sets, modes_to_compute, s::DTAHAS, factory::JuMP.OptimizerFactory,
                         set_variables::AbstractVector{<:SetProg.AbstractVariable} = map(cv->Ellipsoid(InteriorPoint(cv[1])),
                                                                                         chebyshevcenter.(stateset.(s.modes)));
                         y=_vars(s),
                         cone=SOSCone(),
                         λ=Dict{transitiontype(s), Float64}(),
                         enabled = 1:nstates(s),
                         volume_heuristic = nth_root,
                         verbose=1)
    model = SOSModel(factory)
    #set_vrefs = @variable(model, [q in modes_to_compute], Ellipsoid(point=h[q]))
    # FIXME calls JuMP.variable_type(model, Ellipsoid(point=h[q])) then
    # complains that q is not defined
    set_vrefs = Dict(state => @variable(model, variable_type=set_variables[state]) for state in modes_to_compute)

    @objective(model, Max,
               sum(set -> volume_heuristic(volume(set)), values(set_vrefs)))

    λouts = Dict{transitiontype(s), JuMP.AffExpr}()

    for q in modes_to_compute
        @constraint(model, set_vrefs[q] ⊆ stateset(s, q))
        # Invariance constraint
        for t in out_transitions(s, q)
            λin = get(λ, t, nothing)
            target_state = target(s, t)
            if target_state in enabled
                source_set = set_vrefs[q]
                target_set = target_state in modes_to_compute ? set_vrefs[target_state] : sets[target_state]
                r = s.resetmaps[symbol(s, t)]
                λouts[t] = 1.0
                @constraint(model, r.A * source_set ⊆ r.E * target_set,
                            S_procedure_scaling = λin)
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

    for q in modes_to_compute
        sets[q] = JuMP.value(set_vrefs[q])
    end

    return sets
end

function invariant_sets(s::DTAHAS, args...; kws...)
    nmodes = nstates(s)
    sets = Vector{SetProg.Sets.AbstractSet{Float64}}(undef, nmodes)
    return invariant_sets!(sets, 1:nmodes, s, args...; kws...)
end
function invariant_sets!(sets, modes_to_compute, s::DTAHCS, args...; kws...)
    return invariant_sets!(sets, modes_to_compute, algebraiclift(s), args...;
                           kws...)
end
function invariant_sets(s::DTAHCS, args...; kws...)
    return invariant_sets(algebraiclift(s), args...; kws...)
end

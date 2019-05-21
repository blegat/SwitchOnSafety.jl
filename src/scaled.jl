struct ScaledHybridSystem{T, H <: Union{DiscreteSwitchedLinearSystem, ConstrainedDiscreteSwitchedLinearSystem}} <: HybridSystems.AbstractHybridSystem
    system::H
    γ::T
end

# The hybrid system both acts like an automaton and a system

# Automaton
for f in [:state_property, :transition_property, :states, :nstates, :transitiontype, :transitions, :ntransitions, :source, :event, :target, :has_transition, :rem_transition!, :in_transitions, :out_transitions]
    @eval begin
        HybridSystems.$f(hs::ScaledHybridSystem, args...) = HybridSystems.$f(hs.system, args...)
    end
end
for f in [:state_property_type, :transition_property_type]
    @eval begin
        HybridSystems.$f(::Type{<:ScaledHybridSystem{T, H}}, args...) where {T, H} = HybridSystems.$f(H, args...)
    end
end

# System
MathematicalSystems.statedim(hs::ScaledHybridSystem, u) = MathematicalSystems.statedim(hs.system, u)
MathematicalSystems.stateset(hs::ScaledHybridSystem, u) = MathematicalSystems.stateset(hs.system, u)
function scale(lm::LinearMap, γ)
    return LinearMap(lm.A / γ)
end
function scale(s::LinearDiscreteSystem, γ)
    return LinearDiscreteSystem(s.A / γ)
end


function HybridSystems.resetmap(hs::ScaledHybridSystem, t)
    scale(HybridSystems.resetmap(hs.system, t), hs.γ)
end

# TODO should be able to remove it with SetProg
MultivariatePolynomials.variables(s::ScaledHybridSystem, state::Int) = MultivariatePolynomials.variables(s.system, state)

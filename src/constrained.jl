export ConstrainedDiscretePeriodicSwitching, ConstrainedDiscreteSwitchedSystem

type ConstrainedDiscretePeriodicSwitching <: AbstractPeriodicSwitching
    s::AbstractDiscreteSwitchedSystem # Cannot say DiscreteSwitchedSystem as it would be a circular type declaration https://github.com/JuliaLang/julia/issues/269
    period::Vector{Edge}
    growthrate::Float64
end

type ConstrainedDiscreteSwitchedSystem <: AbstractDiscreteSwitchedSystem
    A::Vector
    # Automaton with constrains
    G::DiGraph
    # matrix corresponding to edge is A[σ[edge]]
    σ::Dict{Edge, Int}
    eid::Dict{Edge, Int}
    # Dimension of the nodes : n[i] is the dimension of the `i`th node
    n::Vector{Int}
    x::Vector{Vector{PolyVar{true}}}
    lb::Float64
    ub::Float64
    # There will typically only be lyapunov for small d so a dictionary would be overkill
    lyaps::Vector{Nullable{Lyapunov}}
    smp::Nullable{ConstrainedDiscretePeriodicSwitching}
    function ConstrainedDiscreteSwitchedSystem(A::Vector, G::DiGraph, σ::Dict{Edge, Int})
        @assert !isempty(A) "Needs at least one matrix in the system"
        @assert length(σ) == ne(G) "Number of labels different that number of edges"
        @assert 1 <= maximum(values(σ)) <= length(A) "Invalid labels for the edges"
        n = zeros(Int, nv(G))
        eid = Dict{Edge, Int}()
        neid = 0
        for e in edges(G)
            neid += 1
            eid[e] = neid
            M = A[σ[e]]
            @assert isa(M, AbstractMatrix) "One of the matrices is of invalid type: $(typeof(M))"
            for (u, nu) in ((e.first, size(M, 2)), (e.second, size(M, 1)))
                if n[u] == 0
                    n[u] = nu
                else
                    @assert n[u] == nu "Inconsistent matrix dimensions at node $u"
                end
            end
        end
        y = Vector{Vector{PolyVar{true}}}(nv(G))
        for v in 1:nv(G)
            @polyvar x[1:n[v]]
            y[v] = x
        end
        new(A, G, σ, eid, n, y, 0, Inf, Nullable{Vector{Lyapunov}}[], nothing)
    end
end

nnodes(s::ConstrainedDiscreteSwitchedSystem) = nv(s.G)
MultivariatePolynomials.vars(s::ConstrainedDiscreteSwitchedSystem, i::Int) = s.x[i]

ρA(s::ConstrainedDiscreteSwitchedSystem) = ρ(adjacency_matrix(s.G))

startnode(edge::Edge) = edge.first
endnode(edge::Edge) = edge.second

dynamicfor(s::AbstractSwitchedSystem, edge::Edge) = dynamicfor(s, s.σ[edge])

function state(s::ConstrainedDiscreteSwitchedSystem, edge::Edge, forward=true)
    forward ? edge.second : edge.first
end

function modes(s::ConstrainedDiscreteSwitchedSystem, v::Int, forward=true)
    if forward
        out_edges(s.G, v)
    else
        in_edges(s.G, v)
    end
end

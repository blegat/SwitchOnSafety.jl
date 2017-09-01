export ConstrainedDiscretePeriodicSwitching, ConstrainedDiscreteSwitchedSystem

struct ConstrainedDiscretePeriodicSwitching <: AbstractPeriodicSwitching
    s::AbstractDiscreteSwitchedSystem # Cannot say DiscreteSwitchedSystem as it would be a circular type declaration https://github.com/JuliaLang/julia/issues/269
    period::Vector{Edge}
    growthrate::Float64
end

mutable struct ConstrainedDiscreteSwitchedSystem <: AbstractDiscreteSwitchedSystem
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
            for (u, nu) in ((e.src, size(M, 2)), (e.dst, size(M, 1)))
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
MultivariatePolynomials.variables(s::ConstrainedDiscreteSwitchedSystem, i::Int) = s.x[i]

ρA(s::ConstrainedDiscreteSwitchedSystem) = ρ(adjacency_matrix(s.G))

startnode(edge::Edge) = edge.src
endnode(edge::Edge) = edge.dst

dynamicfor(s::AbstractSwitchedSystem, edge::Edge) = dynamicfor(s, s.σ[edge])

function state(s::ConstrainedDiscreteSwitchedSystem, edge::Edge, forward=true)
    forward ? edge.dst : edge.src
end

function modes(s::ConstrainedDiscreteSwitchedSystem, v::Int, forward=true)
    if forward
        Edge.(v, out_neighbors(s.G, v))
    else
        Edge.(in_neighbors(s.G, v), v)
    end
end

function nullsmp(s::ConstrainedDiscreteSwitchedSystem)
    Nullable{ConstrainedDiscretePeriodicSwitching}()
end
function buildsmp(s::ConstrainedDiscreteSwitchedSystem, seq, growthrate, dt)
    Nullable{ConstrainedDiscretePeriodicSwitching}(ConstrainedDiscretePeriodicSwitching(s, seq, growthrate))
end

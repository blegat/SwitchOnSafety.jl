using LightGraphs
export sosextractsequence

function findapprox(a, x, atol=1e-4)
    for i in eachindex(a)
        if isapprox(a[i][1], x, atol=1e-4)
            return i
        end
    end
    return 0
end

function pushapprox!(a, x, id, args...)
    i = findapprox(a, x, args...)
    if iszero(i)
        push!(a, (x, id))
        id
    else
        a[i][2]
    end
end

function startstates(s::AbstractDiscreteSwitchedSystem, u, G, B,lyap, ranktol=1e-4)
    a = Tuple{Vector{Float64}, Int}[]
    for edge in modes(s, u)
        σ = s.σ[edge]
        μ = lyap.dual[σ]
        X = monomials(variables(μ), d)
        ν = matmeasure(μ, X)
        for x in extractatoms(ν, ranktol).support
            x /= norm(x)
            i = pushapprox!(a, x, length(G)+1)
            if i == length(G)+1
                push!(G, [(σ, endnode(edge))])
                push!(B, x)
            else
                push!(G[i], (σ, endnode(edge)))
            end
        end
    end
    a
end

# Atom extraction
function sosextractsequence(s::AbstractDiscreteSwitchedSystem, d::Integer; solver::AbstractMathProgSolver=JuMP.UnsetSolver(), tol=1e-5)::Nullable{periodicswitchingtype(s)}
    lyap = getlyap(s, d; solver=solver, tol=tol)
    G = Vector{Tuple{Int, Int}}[] # G[u] = list of edges (σ, v) going out of u
    B = Vector{Float64}[] # B[u] = values of the state at node u
    GB = map(u -> startstates(s, u, G, B, lyap), 1:nnodes(s))
    g = DiGraph(length(G))
    w = Matrix{Float64}(nv(g), nv(g))
    ij2σ = Dict{NTuple{2, Int}, Int}()
    for u in eachindex(G)
        x = B[u]
        for (i, (σ, V)) in enumerate(G[u])
            y = s.A[σ] * x
            λ = norm(y)
            y /= λ
            j = findapprox(GB[V], y)
            if j == 0
                return nothing
            end
            v = GB[V][j][2]
            add_edge!(g, v)
            w[u, v] = λ
            ij2σ[(i, j)] = σ
        end
    end
    c, λmin = karp_minimum_cycle_mean(g, w)
    periodicswitching(s, c)
    @show c
    @show λmin
    period = map(ij2σ[(c[i], c[(i % length(c)) + 1])], eachindex(c))

    periodicswitching()
end

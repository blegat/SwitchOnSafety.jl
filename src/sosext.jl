using LightGraphs
export sosextractcycle

function dist(X, Y)
    maximum(x -> abs(x[1] - x[2]) / max(x[1], x[2]), zip(X, Y))
end

function findapprox(a, x, tol)
    _a(i) = dot(a[i][1], x)
    for i in eachindex(a)
        #if dist # isapprox(a[i][1], x, rtol=tol) # rtol is too restrictive since [0.999, 0.001] should be the same as [0.999, 0.002] but rtol says no
        if isapprox(_a(i), 1., rtol=tol)
            return i
        end
    end
    return 0
end

function pushapprox!(a, x, id, tol)
    i = findapprox(a, x, tol)
    if iszero(i)
        push!(a, (x, id))
        id
    else
        a[i][2]
    end
end

σfor(s, edge::Int) = edge
σfor(s, edge::Edge) = s.σ[edge]

function extractstates(s::AbstractDiscreteSwitchedSystem, d, edge, dual, ranktol)
    σ = σfor(s, edge)
    measurefor(dual, s, edge)
    μ = dual[σ]
    X = monomials(variables(μ), d)
    ν = matmeasure(μ, X)
    atoms = extractatoms(ν, ranktol)
    if isnull(atoms)
        Vector{Float64}[]
    else
        get(atoms).support
    end
end

function startstates(s::AbstractDiscreteSwitchedSystem, edgestates, G, B, disttol)
    a = Tuple{Vector{Float64}, Int}[]
    for (edge, states) in edgestates
        for x in states
            x /= norm(x)
            i = pushapprox!(a, x, length(G)+1, disttol)
            σ = σfor(s, edge)
            if i == length(G)+1
                push!(G, [(σ, edge)])
                push!(B, x)
            else
                push!(G[i], (σ, edge))
            end
        end
    end
    a
end

# Atom extraction
function sosextractcycle(s::AbstractDiscreteSwitchedSystem, dual, d::Integer; ranktols=1e-5, disttols=1e-5)
    smp = Nullable{periodicswitchingtype(s)}(nothing)
    for ranktol in ranktols
        # This part is the more costly since it does atom extraction
        # It is run only once for each disttols which is nice
        edgestates = map(u -> map(edge -> (edge, extractstates(s, d, edge, dual, ranktol)), modes(s, u)), 1:nnodes(s))

        for disttol in disttols
            G = Vector{Tuple{Int, edgetype(s)}}[] # G[u] = list of edges (σ, v) going out of u
            B = Vector{Float64}[] # B[u] = values of the state at node u
            GB = map(u -> startstates(s, edgestates[u], G, B, disttol), 1:nnodes(s))

            # Try grouping with different distances

            g = DiGraph(length(G))
            w = fill(Inf, nv(g), nv(g))
            ij2edge = Dict{NTuple{2, Int}, edgetype(s)}()
            for u in eachindex(G)
                x = B[u]
                for (σ, edge) in G[u]
                    y = s.A[σ] * x
                    λ = norm(y)
                    y /= λ
                    j = findapprox(GB[endnode(edge)], y, disttol)
                    if !iszero(j)
                        v = GB[endnode(edge)][j][2]
                        add_edge!(g, u, v)
                        # we use a minimum cycle arithmetic mean algo but we need the maximum cycle geometric mean
                        w[u, v] = -log(λ)
                        ij2edge[(u, v)] = edge
                    end
                end
            end

            c, λmin = karp_minimum_cycle_mean(g, w)

            if !isempty(c)
                period = map(i -> ij2edge[(c[i], c[(i % length(c)) + 1])], eachindex(c))

                newsmp = periodicswitching(s, period)
                if isnull(smp) || isbetter(newsmp, get(smp))
                    smp = Nullable{periodicswitchingtype(s)}(newsmp)
                end
            end
        end
    end

    if !isnull(smp)
        updatesmp!(s, get(smp))
    end

    smp
end
function sosextractcycle(s::AbstractDiscreteSwitchedSystem, d::Integer; solver::AbstractMathProgSolver=JuMP.UnsetSolver(), tol=1e-5, ranktols=tol, disttols=tol)::Nullable{periodicswitchingtype(s)}
    lyap = getlyap(s, d; solver=solver, tol=tol)
    sosextractcycle(s, lyap.dual, d; ranktols=ranktols, disttols=disttols)
end

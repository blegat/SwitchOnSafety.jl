export pradius, pradiusb

function pradiusk(As, p, k, pnorm)
    mean(map(A->norm(A, pnorm)^p, As))^(1/(p*k))
end

function checkeven{A}(p, algo::Type{Val{A}})
    if !iseven(p)
        throw(ArgumentError("Odd p is not supported yet for pradius computation with $A algorithm"))
    end
end

function pradius(s::DiscreteSwitchedSystem, p, lift::Function; pnorm=Inf, ɛ=1e-2, forceub=false)
    ρpmpp = ρ(sum((A -> lift(A, p)).(s.A)))
    if forceub
        ρpmpp^(1/p)
    else
        (ρpmpp / ρA(s))^(1/p)
    end
end

function kozyakinlift(A, edge::Edge)
    kron(sparse([edge.first], [edge.second], [1]), A)
end
function pradius(s::ConstrainedDiscreteSwitchedSystem, p, lift::Function; pnorm=Inf, ɛ=1e-2, forceub=false)
    Alifted = (A -> lift(A, p)).(s.A)
    ρpmpp = ρ(sum(kozyakinlift(Alifted[s.σ[edge]], edge) for edge in edges(s.G)))
    if forceub
        ρpmpp^(1/p)
    else
        (ρpmpp / ρA(s))^(1/p)
    end
end

function pradius(s::AbstractDiscreteSwitchedSystem, p, algo::Type{Val{:VeroneseLift}}; pnorm=Inf, ɛ=1e-2, forceub=false)
    checkeven(p, algo)
    pradius(s, p, veroneselift, pnorm=pnorm, ɛ=ɛ, forceub=forceub)
end
kroneckerlift(A::AbstractMatrix, p::Integer) = kronpow(veroneselift(A, 2), div(p, 2))
function pradius(s::AbstractDiscreteSwitchedSystem, p, algo::Type{Val{:KroneckerLift}}; pnorm=Inf, ɛ=1e-2, forceub=false)
    checkeven(p, algo)
    pradius(s, p, kroneckerlift, pnorm=pnorm, ɛ=ɛ, forceub=forceub)
end

function pradius(s::DiscreteSwitchedSystem, p, algo::Type{Val{:BruteForce}}; pnorm=Inf, ɛ=1e-2, forceub=false)
    As = s.A
    ρp = Float64[]
    k = 0
    while k < 2 || !isapprox(ρp[end-1], ρp[end], rtol=ɛ)
        k += 1
        m = length(s.A)
        l = length(As)
        Asnew = similar(As, l * m)
        for (i,A) in enumerate(s.A)
            for (j,B) in enumerate(As)
                Asnew[(i-1)*l+j] = A*B
            end
        end
        As = Asnew
        push!(ρp, pradiusk(As, p, k, pnorm))
    end
    if forceub
        ρp[end] * length(s.A)^(1/p)
    else
        ρp[end]
    end
end

# The Inf- and 1-norm of a matrices are easier to compute than the 2-norm
function pradius(s::AbstractDiscreteSwitchedSystem, p, algo::Symbol=:VeroneseLift; pnorm=Inf, ɛ=1e-2, forceub=false)
    pradius(s, p, Val{algo}, pnorm=pnorm, ɛ=ɛ, forceub=forceub)
end

function pradiusb(s::DiscreteSwitchedSystem, p, algo::Symbol=:VeroneseLift)
    @assert algo in [:VeroneseLift, :KroneckerLift] "p-radius algo needs to be exact to compute bounds on JSR"
    ρpmp = pradius(s, p, algo, forceub=true)
    ρp = ρpmp / ρA(s)^(1/p)
    updateb!(s, ρp, ρpmp)
end

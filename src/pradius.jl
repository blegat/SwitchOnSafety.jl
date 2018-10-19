export pradius, pradiusb

function pradiusk(As, p, k, pnorm)
    mean(map(A->opnorm(A, pnorm)^p, As))^(1/(p*k))
end

function pradius(s::DiscreteSwitchedLinearSystem, p, lift::Function; pnorm=Inf, ɛ=1e-2, forceub=false)
    ρpmpp = ρ(sum((rm -> lift(rm.A, p)).(s.resetmaps)))
    if forceub
        ρpmpp^(1/p)
    else
        (ρpmpp / ρA(s))^(1/p)
    end
end

function kozyakinlift(s, t, As)
    n = nstates(s)
    A = As[symbol(s, t)]
    kron(sparse([source(s, t)], [target(s, t)], [one(eltype(A))], n, n), A)
end
function pradius(s::ConstrainedDiscreteSwitchedLinearSystem, p, lift::Function; pnorm=Inf, ɛ=1e-2, forceub=false)
    Alifted = (rm -> lift(rm.A, p)).(s.resetmaps)
    ρpmpp = ρ(sum(kozyakinlift(s, t, Alifted) for t in transitions(s)))
    if forceub
        ρpmpp^(1/p)
    else
        (ρpmpp / ρA(s))^(1/p)
    end
end

abstract type PRadiusAlgorithm end
abstract type ExactPRadiusAlgorithm <: PRadiusAlgorithm end
function checkeven(p, algo::ExactPRadiusAlgorithm)
    if !iseven(p)
        throw(ArgumentError("Odd p is not supported yet for pradius computation with $algo algorithm"))
    end
end

export VeroneseLift, KroneckerLift, BruteForce
struct VeroneseLift <: ExactPRadiusAlgorithm end
struct KroneckerLift <: ExactPRadiusAlgorithm end
struct BruteForce <: PRadiusAlgorithm end

liftfunction(::VeroneseLift) = veroneselift
liftfunction(::KroneckerLift) = kroneckerlift
function kroneckerlift(A::AbstractMatrix, p::Integer)
    @assert iseven(p)
    kronpow(veroneselift(A, 2), div(p, 2))
end


function pradius(s::AbstractDiscreteSwitchedSystem, p, algo::ExactPRadiusAlgorithm=VeroneseLift(); kws...)
    checkeven(p, algo)
    pradius(s, p, liftfunction(algo); kws...)
end
# The Inf- and 1-norm of a matrices are easier to compute than the 2-norm
function pradius(s::DiscreteSwitchedLinearSystem, p, ::BruteForce; pnorm=Inf, ɛ=1e-2, forceub=false)
    As = map(rm -> rm.A, s.resetmaps)
    Ascur = As
    ρp = Float64[]
    m = ntransitions(s)
    k = 0
    while k < 2 || !isapprox(ρp[end-1], ρp[end], rtol=ɛ)
        k += 1
        l = length(Ascur)
        Asnew = similar(Ascur, l * m)
        for (i,A) in enumerate(As)
            for (j,B) in enumerate(Ascur)
                Asnew[(i-1)*l+j] = A*B
            end
        end
        Ascur = Asnew
        push!(ρp, pradiusk(Ascur, p, k, pnorm))
    end
    if forceub
        ρp[end] * length(As)^(1/p)
    else
        ρp[end]
    end
end

# p-radius algo needs to be exact to compute bounds on JSR
function pradiusb(s::AbstractDiscreteSwitchedSystem, p, algo::ExactPRadiusAlgorithm=VeroneseLift())
    ρpmp = pradius(s, p, algo, forceub=true)
    ρp = ρpmp / ρA(s)^(1/p)
    updateb!(s, ρp, ρpmp)
end

using MultivariatePolynomials
using JuMP

using RecipesBase

export Ellipsoid

struct Ellipsoid{T}
    Q::Matrix{T}
    c::Vector{T}
end

@recipe function f(ell::Ellipsoid)
    @assert Base.LinAlg.checksquare(ell.Q) == 2
    αs = linspace(0, 2π, 1024)
    ps = [[cos(α), sin(α)] for α in αs]
    r = [sqrt(dot(p, ell.Q * p)) for p in ps]
    seriestype --> :shape
    legend --> false
    ell.c[1] .+ cos.(αs) ./ r, ell.c[2] .+ sin.(αs) ./ r
end

struct LiftedEllipsoid{T}
    P::Matrix{T}
end

function LiftedEllipsoid(ell::Ellipsoid)
    md = ell.Q*ell.c
    δ = ell.c'*md-1
    d = -md
    D = ell.Q
    P = [δ d'
         d D]
    LiftedEllipsoid(P)
end

Base.convert(::Type{Ellipsoid{T}}, ell::LiftedEllipsoid) where T = convert(Ellipsoid{T}, Ellipsoid(ell))
function Ellipsoid(ell::LiftedEllipsoid)
    # P is
    # λ * [c'Qc-1  -c'Q
    #         -Qc   Q]
    # Let P be [β b'; b B]
    # We have
    # β = λ c'Qc - λ
    # b = -λ Qc <=> Q^{-1/2}b = -λ Q^{1/2} c
    # hence
    # λ c'Qc = β + λ
    # λ^2 c'Qc = b'Q^{-1}b = λ b'B^{-1}b <=> λ c'Qc = b'B^{-1}b
    # Hence λ = b'B^{-1}b - β
    n = LinAlg.checksquare(ell.P)-1
    ix = 1+(1:n)
    β = ell.P[1, 1]
    b = ell.P[1, ix]
    B = ell.P[ix, ix]
    λ = dot(b, B \ b) - β
    c = -(B \ b)
    Q = B / λ
    Ellipsoid(Q, c)
end

@recipe function f(ell::LiftedEllipsoid)
    Ellipsoid(ell)
end

struct ConeLyap{T, P<:AbstractPolynomial{T}, S}
    p::P
    Q::Matrix{S}
    b::Vector{S}
    β::S
    c::Vector{Float64}
    H::Matrix{Float64}
    vol::S
    #L::Matrix{JuMP.Variable}
    #λinv::Vector{JuMP.Variable}
end

ConeLyap(p::P, Q::Matrix{S}, b::Vector{S}, β::S, c, H, vol::S) where {T, P<:AbstractPolynomial{T}, S} = ConeLyap{T, P, S}(p, Q, b, β, c, H, vol)
JuMP.resultvalue(p::ConeLyap) = ConeLyap(JuMP.resultvalue(p.p), JuMP.resultvalue(p.Q), JuMP.resultvalue(p.b), getvalue(p.β), p.c, p.H, JuMP.resultvalue(p.vol))

ellipsoid(p::ConeLyap{T, P, JuMP.Variable}) where {T, P<:AbstractPolynomial{T}} = ellipsoid(JuMP.resultvalue(p))

function _HPH(D, d, δ, H)
    P = [δ d'
         d D]
    HPH = H * P * H
end
function ellipsoid(p::ConeLyap)
    LiftedEllipsoid(inv(_HPH(p.Q, p.b, p.β, p.H)))
end

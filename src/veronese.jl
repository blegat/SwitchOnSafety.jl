export permanent, veroneselift

using Combinatorics

function permanent(A, perms)
    n = LinearAlgebra.checksquare(A)
    sum(σ -> prod(i -> A[i,σ[i]], 1:n), perms)
end

function μ(combi)
    prod = 1
    last = 0
    cur = 1
    for i in combi
        if i == last
            cur += 1
        else
            prod *= factorial(cur)
            last = i
            cur = 1
        end
    end
    prod *= factorial(cur)
    return prod
end

"""
   veroneselift_explicit(A::AbstractMatrix, d::Integer)

Computes the Veronese lift of degree `d` of `A` using the explicit formula
```math
(A^{[d]})_{\\alpha\\beta} = \\frac{\\text{per}\\,A(\\alpha, \\beta)}{\\sqrt{\\mu(\\alpha)\\mu(\\beta)}}
```
given in equation (8) of [PJ08]. Note that this is a lot slower than
[`veroneselift`](@ref) and is implemented mainly for educational purpose.

* [PJ08] P. Parrilo and A. Jadbabaie.
*Approximation of the joint spectral radius using sum of squares*.
Linear Algebra and its Applications, Elsevier, **2008**, 428, 2385-2402
"""
function veroneselift_explicit(A::AbstractMatrix, d::Integer)
    d > 0 || throw(ArgumentError("`d` should be a positive integer, got $d"))

    n = LinearAlgebra.checksquare(A)
    N = binomial(n+d-1, d)

    Ad = Matrix{Float64}(undef, N, N)

    # Computing permutations is the most expensive part so we precompute it
    perms = collect(permutations(1:d))
    for (i,α) in enumerate(with_replacement_combinations(1:n, d))
        for (j,β) in enumerate(with_replacement_combinations(1:n, d))
            Ad[i, j] = permanent(A[α, β], perms) / sqrt(μ(α) * μ(β))
        end
    end

    return Ad
end


function transform_combination(α, m)
    β = zeros(Int, m)
    for i in α
        β[i] += 1
    end
    return β
end

"""
   veroneselift(A::AbstractMatrix, d::Integer)

Computes the Veronese lift of degree `d` of `A` using `DynamicPolynomials`.
"""
function veroneselift(A::AbstractMatrix, d::Integer)
    d > 0 || throw(ArgumentError("`d` should be a positive integer, got $d"))
    m, n = size(A)
    @polyvar x[1:n]
    Ax = A * x
    M = binomial(m + d - 1, d)
    N = binomial(n + d - 1, d)
    Ad = Matrix{Float64}(undef, M, N)
    df = factorial(d)
    scaling(m) = sqrt(df / prod(factorial, exponents(m)))
    X = monomials(x, d)
    # Scaling for the scaled monomial basis
    col_scaling = Float64[scaling(m) for m in X]
    if m == n
        Y = X
        row_scaling = col_scaling
    else
        # Only exponents matter
        @polyvar y[1:m]
        Y = monomials(y, d)
        row_scaling = Float64[scaling(m) for m in Y]
    end
    for (i, mono) in enumerate(Y)
        α = exponents(mono)
        Axα = prod(i -> Ax[i]^α[i], 1:m)
        Ad[i, :] = coefficients(Axα, X) .* row_scaling[i] ./ col_scaling
    end
    return Ad
end

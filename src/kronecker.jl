# Computes A^(⊗p) = A ⊗ A ⊗ ... ⊗ A
# where ⊗ is the Kronecker product
# Since p > 0, B will be kron'd at least once so it will be a Matrix{T}
# This is not easy to guess for Julia inference so I annotate the return type
function (kronpow(A::MT, p)::MT) where MT <: AbstractMatrix
    @assert p > 0
    B = 1
    C = A
    bitmap = 1
    while bitmap <= p
        if (bitmap & p) != 0
            B = kron(B, C)
        end
        C = kron(C, C)
        bitmap <<= 1
    end
    B
end

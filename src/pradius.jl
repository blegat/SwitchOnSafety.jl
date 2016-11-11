export pradius, pradiusb

function pradiusk(As, p, k, pnorm)
    mean(map(A->norm(A, pnorm)^p, As))^(1/(p*k))
end

# The Inf- and 1-norm of a matrices are easier to compute than the 2-norm
function pradius(s::SwitchedSystem, p, algo=:VeroneseLift; pnorm=Inf, ɛ=1e-2, forceub=false)
    if algo == :VeroneseLift
        if !iseven(p)
            throw(ArgumentError("Odd p is not supported yet for pradius computation"))
            end
            ρpmpp = ρ(sum(map(A->veroneselift(A, p), s.A)))
            if forceub
                ρpmpp^(1/p)
            else
                (ρpmpp / length(s.A))^(1/p)
            end
        elseif algo == :BruteForce
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
        else
            throw(ArgumentError("Algorithm $algo not recognized"))
        end
    end

    function pradiusb(s::SwitchedSystem, p)
        ρpmp = pradius(s, p, :VeroneseLift, forceub=true)
        ρp = ρpmp / length(s.A)^(1/p)
        updateb!(s, ρp, ρpmp)
    end

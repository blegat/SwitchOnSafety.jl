function quickb(s::AbstractDiscreteSwitchedSystem, k::Integer=1, clb=true, cub=true)
    lb = 0.
    ub1 = 0.
    ub2 = 0.
    ub∞ = 0.
    for st in states(s)
        for sw in switchings(s, k, st)
            if clb && isperiodic(sw)
                psw = periodicswitching(sw)
                lb = max(lb, psw.growthrate)
                notifyperiodic!(s, psw)
                updatesmp!(s, psw)
            end
            if cub
                ub1 = max(ub1, opnorm(sw.A, 1)^(1/k))
                ub2 = max(ub2, opnorm(sw.A, 2)^(1/k))
                ub∞ = max(ub∞, opnorm(sw.A, Inf)^(1/k))
            end
        end
    end
    ub = min(ub1, ub2, ub∞)
    cub && updateub!(s, ub)
    lb, ub
end

function quickb(s::AbstractDiscreteSwitchedSystem, K::AbstractVector, args...)
    reduce((x, y) -> (max(x[1],y[1]), min(x[2], y[2])), quickb.(s, K, args...), init = (0., Inf))
end

ρA(s::DiscreteSwitchedLinearSystem) = ntransitions(s)

function quicklb(s::DiscreteSwitchedLinearSystem, k=1)
    quickb(s, k, true, false)[1]
end

function quickub(s::AbstractDiscreteSwitchedSystem, k=1)
    quickb(s, k, false, true)[2]
end

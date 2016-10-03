# Similar to JuMP/test/solvers.jl

function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

mos = try_import(:Mosek)
scs = false # try_import(:SCS) # It does not work

# Semidefinite solvers
sdp_solvers = Any[]
mos && push!(sdp_solvers, Mosek.MosekSolver(LOG=0))
scs && push!(sdp_solvers, SCS.SCSSolver(verbose=0))

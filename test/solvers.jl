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
ismosek(solver) = contains(string(typeof(solver)),"MosekSolver")
csd = try_import(:CSDP)
iscsdp(solver) = contains(string(typeof(solver)),"CSDP")
scs = false && try_import(:SCS) # It does not work
isscs(solver) = contains(string(typeof(solver)),"SCSSolver")
ipt = try_import(:Ipopt)

# Semidefinite solvers
sdp_solvers = Any[]
mos && push!(sdp_solvers, Mosek.MosekSolver(LOG=0))
csd && push!(sdp_solvers, CSDP.CSDPSolver(printlevel=0))
scs && push!(sdp_solvers, SCS.SCSSolver(verbose=0))

# Bilinear LP solvers
blp_solvers = Any[]
ipt && push!(blp_solvers, Ipopt.IpoptSolver(print_level=0))

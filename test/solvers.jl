# Similar to JuMP/test/solvers.jl

function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

mos = false && try_import(:Mosek)
if mos
    mossolver = () -> Mosek.MosekInstance(LOG=0)
else
    mossolver = () -> nothing
end
ismosek(solver) = solver === mossolver
csd = false && try_import(:CSDP)
if csd
    csdsolver = () -> CSDP.CSDPInstance(printlevel=0)
else
    csdsolver = () -> nothing
end
iscsdp(solver) = solver === csdsolver
sda = try_import(:SDPA)
if sda
    sdasolver = () -> SDPA.SDPAInstance()
else
    sdasolver = () -> nothing
end
issdpa(solver) = solver === sdasolver
scs = false && try_import(:SCS) # It does not work
isscs(solver) = false
ipt = false && try_import(:Ipopt)

# Semidefinite solvers
sdp_solvers = Any[]
mos && push!(sdp_solvers, mossolver)
csd && push!(sdp_solvers, csdsolver)
sda && push!(sdp_solvers, sdasolver)
scs && push!(sdp_solvers, scssolver)

# Bilinear LP solvers
blp_solvers = Any[]
ipt && push!(blp_solvers, Ipopt.IpoptSolver(print_level=0))

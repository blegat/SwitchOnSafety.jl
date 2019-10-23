# Similar to JuMP/test/solvers.jl

using JuMP

function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

mos = try_import(:MosekTools)
if mos
    mos_factory = with_optimizer(MosekTools.Mosek.Optimizer, QUIET=true)
else
    mos_factory = nothing
end
ismosek(factory) = factory === mos_factory
csd = try_import(:CSDP)
if csd
    csd_factory = with_optimizer(CSDP.Optimizer, printlevel=0)
else
    csd_factory = nothing
end
iscsdp(solver) = solver === csd_factory
sda = try_import(:SDPA)
if sda
    sda_factory = with_optimizer(SDPA.Optimizer)
else
    sda_factory = nothing
end
issdpa(solver) = solver === sda_factory
scs = false && try_import(:SCS) # It does not work
isscs(solver) = false
ipt = false && try_import(:Ipopt)

# Semidefinite solvers
sdp_factories = Any[]
mos && push!(sdp_factories, mos_factory)
csd && push!(sdp_factories, csd_factory)
sda && push!(sdp_factories, sda_factory)
#scs && push!(sdp_factories, scs_factory)

# Bilinear LP solvers
blp_factories = Any[]
ipt && push!(blp_factories, with_optimizer(Ipopt.Optimizer, print_level=0))

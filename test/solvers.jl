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
    mos_optimizer_constructor = optimizer_with_attributes(MosekTools.Mosek.Optimizer, "QUIET" => true)
else
    mos_optimizer_constructor = nothing
end
ismosek(optimizer_constructor) = optimizer_constructor === mos_optimizer_constructor
csd = try_import(:CSDP)
if csd
    csd_optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, "printlevel" => 0)
else
    csd_optimizer_constructor = nothing
end
iscsdp(solver) = solver === csd_optimizer_constructor
sda = try_import(:SDPA)
if sda
    sda_optimizer_constructor = SDPA.Optimizer
else
    sda_optimizer_constructor = nothing
end
issdpa(solver) = solver === sda_optimizer_constructor
scs = false && try_import(:SCS) # It does not work
isscs(solver) = false
ipt = false && try_import(:Ipopt)
eco = try_import(:ECOS)
if eco
    eco_optimizer_constructor = optimizer_with_attributes(ECOS.Optimizer, MOI.Silent() => true)
else
    eco_optimizer_constructor = nothing
end


# Accurate SOC solvers
soc_optimizer_constructors = Any[]
eco && push!(soc_optimizer_constructors, eco_optimizer_constructor)

# Semidefinite solvers
sdp_factories = Any[]
mos && push!(sdp_factories, mos_optimizer_constructor)
csd && push!(sdp_factories, csd_optimizer_constructor)
sda && push!(sdp_factories, sda_optimizer_constructor)
#scs && push!(sdp_factories, scs_optimizer_constructor)

# Bilinear LP solvers
blp_factories = Any[]
ipt && push!(blp_factories, optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))

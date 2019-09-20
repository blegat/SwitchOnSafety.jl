@testset "Identity dynamical system in square" begin
    include("identity.jl")
end
@testset "Linear dynamical system in square" begin
    include("lti.jl")
end
@testset "Controlled linear dynamical system in square" begin
    include("controlled.jl")
end
include("cis.jl")

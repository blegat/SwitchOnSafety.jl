@testset "Unconstrained" begin
    include("finiteness_conjecture_counterexample.jl")
    include("double_msd.jl")
    include("AJPR14e54.jl")
    include("AP12e21.jl")
    include("JCG14e61.jl")
    include("JCG14e63.jl")
    include("PJ08e28.jl")
    include("PJ08e54.jl")
end
@testset "Constrained" begin
    include("PEDJ16s4.jl")
end

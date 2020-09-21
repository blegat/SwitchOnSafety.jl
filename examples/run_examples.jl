using Test

const EXAMPLES = ["AJPR14e54.jl"]

@testset "run_examples.jl" begin
    @testset "$(example)" for example in EXAMPLES
        include(example)
    end
end

using Polyhedra
using SwitchOnSafety

@testset "Conitope" begin
    p = polyhedron(Conitope(2, [[1.0, 1.0]]))
    @test nallrays(p) == 0
    @test npoints(p) == 4
end

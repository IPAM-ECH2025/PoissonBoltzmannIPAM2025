using Pkg, JuliaMPBSolver, ExampleJuggler, Test, Markdown

ExampleJuggler.verbose!(true)

notebooks = [
    "ICMPBP-DD-Draft.jl",
]
@testset "Notebooks" begin
    @testscripts(joinpath(@__DIR__, "..", "notebooks"), notebooks)
end

Pkg.test("JuliaMPBSolver")

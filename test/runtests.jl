using Pkg, JuliaMPBSolver, ExampleJuggler, Test, Markdown

ExampleJuggler.verbose!(true)

notebooks = [
    "MPBP-Draft.jl",
    "ICMPBP-Draft.jl",
    "ICMPBP-EndOfHackathon.jl",
    "ICMPBP-DD-Draft.jl",
    "HalfCellAppliedPotential.jl",
    "SymmetricCellSurfaceCharge.jl",
]

@testset "Notebooks" begin
    @testscripts(joinpath(@__DIR__, "..", "notebooks"), notebooks)
end

Pkg.test("JuliaMPBSolver")

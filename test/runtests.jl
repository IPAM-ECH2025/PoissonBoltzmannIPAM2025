using Pkg
using JuliaMPBSolver
using ExampleJuggler
using Test
using Markdown
using JLD2

ExampleJuggler.verbose!(true)

notebooks = [
    "MPBP-Draft.jl",
    "ICMPBP-Draft.jl",
    "ICMPBP-EndOfHackathon.jl",
    "ICMPBP-DD-Draft.jl",
    "HalfCellAppliedPotential.jl",
    "SymmetricCellSurfaceCharge.jl",
]

scripts = [
    "simplecell.jl",
]

# @testset "Notebooks" begin
#     @testscripts(joinpath(@__DIR__, "..", "notebooks"), notebooks)
# end

@testset "Scripts" begin
    @testscripts(joinpath(@__DIR__, "..", "scripts"), scripts)

    for script in scripts
        script_name = script[1:(end - 3)]
        solution, X, nv, ε_r = load(script_name * ".jld2", "solution", "X", "nv", "ε_r")
        solution_reference, X_reference, nv_reference, ε_r_reference = load(script_name * "-reference.jld2", "solution", "X", "nv", "ε_r")

        @test solution ≈ solution_reference
        @test X ≈ X_reference
        @test nv ≈ nv_reference
        @test ε_r ≈ ε_r_reference

        rm(script_name * ".jld2")
    end

end

Pkg.test("JuliaMPBSolver")

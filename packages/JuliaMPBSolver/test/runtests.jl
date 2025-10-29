using JuliaMPBSolver
using Test

@testset "JuliaMPBSolver" begin
  @testset "Units" begin
    @test JuliaMPBSolver.Units.RT(0) == 0
  end
  @testset "Grids" begin
    # TODO: Write unit test here
  end
end
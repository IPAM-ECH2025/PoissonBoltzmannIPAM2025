module PoissonBoltzmannIPAM2025
using DrWatson, PoissonBoltzmannIPAM2025
using PythonPlot
using PythonPlot: pyplot
using LessUnitful
using ExtendableGrids
using VoronoiFVM
using Test
using JuliaMPBSolver.ICMPBP: ICMPBP, SurfaceChargedSymmetricCell, AbstractHalfCell, AbstractSymmetricCell, set_molarity!, calc_cmol, calc_c0mol, calc_χ, get_E, get_φ, get_p, get_c0,
    set_κ!, set_q!, set_φ!, pramp

resultsdir(args...) = projectdir("results", args...)
draftresultsdir(args...) = projectdir("draftresults", args...)

include("plotcells.jl")

export resultsdir, draftresultsdir, plotcells
end

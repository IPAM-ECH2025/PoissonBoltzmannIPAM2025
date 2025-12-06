using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using LinearAlgebra
using Interpolations
using VoronoiFVM
using PythonPlot
using AugmentedPoissonBoltzmann
using JLD2

const nel = 20.0 * AugmentedPoissonBoltzmann.Units.el_surface_density # number of electrons/nm^2 at interfaces
const M_bulk = 1 # (bulk) molarity at center of domain
const E0 = 10AugmentedPoissonBoltzmann.Units.V / AugmentedPoissonBoltzmann.Units.nm # decrement parameter
const a = 5.0 / E0^2 # decrement parameter in χ(E)
const c̄ = 55.508AugmentedPoissonBoltzmann.Units.M # summary molar concentration
const z = [-1, 1]
const c_bulk = [M_bulk / abs(z[1]), M_bulk / abs(z[2])] * AugmentedPoissonBoltzmann.Units.M # bulk  concentrations

# Parameters
user_parameters = AugmentedPoissonBoltzmann.Parameters.UserParameters(
    273.15 + 25 * AugmentedPoissonBoltzmann.Units.K,
    78.49 - 1,
    0.0,
    0.0,
    z,
    c̄,
    c_bulk,
    nel,
    a,
    true,
    false,
    false,
    0.0,
)

computed_parameters =
    AugmentedPoissonBoltzmann.Parameters.ComputedParameters(user_parameters)

# Grid generation
grid_parameters = AugmentedPoissonBoltzmann.Grid.GeometricGrid(
    domain_size = 10.0 * AugmentedPoissonBoltzmann.Units.nm,
    refinement = 4,
    hmin = 1.0e-1 * AugmentedPoissonBoltzmann.Units.nm,
    hmax = 1.0 * AugmentedPoissonBoltzmann.Units.nm,
    use_offset = false,
)

solution, X, nv, ε_r =
    AugmentedPoissonBoltzmann.Equations.create_and_run_full_cell_problem(
    grid_parameters,
    user_parameters,
    computed_parameters,
)

save("simplecell.jld2", "solution", solution, "X", X, "nv", nv, "ε_r", ε_r)

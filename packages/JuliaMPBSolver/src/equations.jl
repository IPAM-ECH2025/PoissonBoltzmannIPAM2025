module Equations

using VoronoiFVM
using ExtendableGrids
using PreallocationTools

using ..Parameters
using ..Units
using ..Grid
using ..Postprocess

function add_dirichlet_bc!(
  system::VoronoiFVM.AbstractSystem,
  field_id::Int,
  face_id::Int,
  value::AbstractFloat,
)
  # Note that because this relies on the equation system being created,
  # method must be called after create_equation_system
  boundary_dirichlet!(system, field_id, face_id, value)
  return nothing
end

function add_boundary_voltage!(
  system::VoronoiFVM.AbstractSystem,
  field_id::Int,
  face_id::Int,
  value::AbstractFloat,
)
  add_dirichlet_bc!(system, field_id, face_id, value)
  return nothing
end

function pin_pressure_value!(
  system::VoronoiFVM.AbstractSystem,
  field_id::Int,
  face_id::Int,
)
  add_dirichlet_bc!(system, field_id, face_id, 0.0)
  return nothing
end

function add_neumann_bc!(
  system::VoronoiFVM.AbstractSystem,
  field_id::Int,
  face_id::Int,
  value::AbstractFloat,
)
  # Note that because this relies on the equation system being created,
  # method must be called after create_equation_system
  boundary_neumann!(system, field_id, face_id, value)
  return nothing
end

function add_boundary_charge!(
  system::VoronoiFVM.AbstractSystem,
  field_id::Int,
  face_id::Int,
  value::AbstractFloat,
)
  add_neumann_bc!(system, field_id, face_id, value)
  return nothing
end

function create_equation_system(
  grid::ExtendableGrid,
  user_parameters::UserParameters,
)
  # Determine the precision of the floating-point numbers from the user parameters
  float_type = eltype(user_parameters.temperature)

  # Create a tmp vector for various function evaluations
  tmp = DiffCache(ones(float_type, length(user_parameters.charge_numbers)))

  function flux!(y, u, edge, data)
    ones = one(eltype(u))
    edge_length = float_type(edgelength(edge))
    electric_field = (u[1, 1] - u[1, 2]) / edge_length

    electric_susceptibility =
      user_parameters.dielectric_susceptibility / sqrt(
        ones +
        user_parameters.electric_susceptibility_decrement_parameter *
        electric_field^2,
      )

    permittivity = (ones + electric_susceptibility) * Units.vacuum_permittivity

    y[1] = permittivity * (u[1, 1] - u[1, 2])

    return nothing
  end

  function reaction!(y, u, node, data)
    local_tmp = get_tmp(tmp, u)
    y[1] = -Postprocess.compute_spacecharge(local_tmp, u[1], user_parameters)
    return nothing
  end

  system = VoronoiFVM.System(
    grid;
    reaction = reaction!,
    flux = flux!,
    species = [1],
    valuetype = float_type,
  )

  return system
end

function solve_equation_system(
  system::VoronoiFVM.AbstractSystem;
  inival = 0.1,
  verbose = "n",
  damp_initial = 0.1,
  maxiters = 1000,
)
  return solve(
    system;
    inival = inival,
    verbose = verbose,
    damp_initial = damp_initial,
    maxiters = maxiters,
  )
end

function create_and_run_full_cell_problem(
  grid_paramters::GeometricGrid,
  user_parameters::UserParameters,
)
  grid = create_full_cell(grid_paramteters)
  system = create_equation_system(grid, user_parameters)

  add_boundary_charge!(system, 1, 2, -Units.elementary_charge)
  add_boundary_charge!(system, 1, 1, Units.elementary_charge)

  solution = solve_equation_system(system)

  return solution
end

export add_boundary_voltage!,
  pin_pressure_value!, add_boundary_charge!, create_equation_system

end

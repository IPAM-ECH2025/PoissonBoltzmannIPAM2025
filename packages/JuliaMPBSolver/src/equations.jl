module Equations

using VoronoiFVM

using ..Parameters

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

function create_equation_system(user_parameters::UserParameters)
  # Determine the precision of the floating-point numbers from the user parameters
  float_type = eltype(user_parameters.temperature)

  # Create a tmp vector for various function evaluations
  tmp = DiffCache(ones(float_type, length(user_parameters.charge_numbers)))

  return nothing
end

export add_boundary_voltage!,
  pin_pressure_value!, add_boundary_charge!, create_equation_system

end

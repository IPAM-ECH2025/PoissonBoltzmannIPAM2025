module Postprocess

using ..Parameters
using ..Units

function compute_mole_fractions!(
  mole_fractions::AbstractVector{<:AbstractFloat},
  electric_potential::AbstractVector{<:AbstractFloat},
  user_parameters::UserParameters,
)
  n_species = length(user_parameters.charge_numbers)

  @assert length(mole_fractions) == n_species "size mismatch between mole_fractions and charge_numbers"

  for i in 1:n_species
    mole_fractions[i] = exp(
      -user_parameters.charge_numbers[i] * electric_potential * F /
      RT(user_parameters.temperature),
    )
  end
end

function compute_concentrations() end

end
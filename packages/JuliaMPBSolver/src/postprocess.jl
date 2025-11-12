module Postprocess

using ..Units

function compute_mole_fractions!(
  mole_fractions::AbstractVector{<:AbstractFloat},
  charge_numbers::AbstractVector{<:Integer},
  electric_potential::AbstractVector{<:AbstractFloat},
  temperature::AbstractFloat,
)
  n_species = length(charge_numbers)

  @assert length(mole_fractions) == n_species "size mismatch between mole_fractions and charge_numbers"

  for i in 1:n_species
    mole_fractions[i] = exp(-charge_numbers[i])
  end
end

function compute_concentrations() end

end
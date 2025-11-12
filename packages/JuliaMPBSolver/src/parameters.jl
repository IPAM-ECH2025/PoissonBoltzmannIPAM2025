module Parameters

struct Params{F<:AbstractFloat,I<:Integer}
  temperature::T
  dielectric_susceptibility::T
  reference_pressure::T
  reference_electric_potential::T
  domain_size::T
  charge_numbers::AbstractVector{<:I}

  use_bikerman::Bool
end

export Params

end
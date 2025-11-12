module Parameters

struct UserParameters{F<:AbstractFloat,I<:Integer}
  temperature::F
  dielectric_susceptibility::F
  reference_pressure::F
  reference_electric_potential::F
  charge_numbers::AbstractVector{<:I}

  use_bikerman_model::Bool
end

export UserParameters

end
using ExtendableGrids
using LessUnitful
using VoronoiFVM

"""
    prepare_grid(n, domain)

Prepare eqidistant simulation grid. Add  interface region in the mid of the
domain in order to fix pressure and provide location for ion conservation
constraints.
"""

function prepare_grid(n, domain)
    if iseven(n)
        n = n + 1
    end
    X = range(domain..., length = n)
    grid = ExtendableGrids.simplexgrid(X)
    xmid = (domain[1] + domain[2]) / 2
    bfacemask!(grid, [xmid], [xmid], 3, tol = 1.0e-10 * ufac"nm")
    return X, grid
end

"""
    mpbpsolve(; n=21, 
                chargenumbers=[-1, 1], 
                domain=[0, 1.0e-10], 
                surfacecharge=[0.16, -0.16], 
                bulkmolarity=0.1)

Solve the modified Poisson-Boltzmann problem with specified bulk molarity.

This function solves the modified Poisson-Boltzmann equation for a 1D domain with specified
surface charges and bulk ion concentrations. The solver uses a finite volume method
with ExtendableGrids and VoronoiFVM.

# Arguments
- `n::Int=21`: Number of grid points along the domain. Default is 21.
- `chargenumbers::Vector{Int}=[-1, 1]`: Charge numbers of the ionic species (e.g., [-1, 1] for monovalent anions and cations).
- `domain::Vector{Float64}=[0, 1.0e-10]`: Domain boundaries in meters [start, end]. Default is [0, 1.0e-10] (0 to 0.1 nm).
- `surfacecharge::Vector{Float64}=[0.16, -0.16]`: Surface charge densities at the domain boundaries in C/m².
- `bulkmolarity::Float64=0.1`: Bulk molarity of the electrolyte solution in mol/L.

# Returns
- `X`: Grid coordinates in m (Vector{Float64})
- `c0`: Solvent concentration in mol/L  (Vector{Float64})
- `c1`: Concentration of first ionic species along the grid  in mol/L (Vector{Float64})
- `c2`: Concentration of second ionic species along the grid  in mol/L(Vector{Float64})

# Examples
```julia
# Solve with default parameters
X, c0, c_anion, c_cation = mpbpsolve()

# Solve with custom parameters
X, c0, c_anion, c_cation = mpbpsolve(
    n=51,
    chargenumbers=[-1, 1],
    domain=[0, 2.0e-10],
    surfacecharge=[0.2, -0.2],
    bulkmolarity=0.15
)
```

# Notes
- The function automatically creates a uniform grid over the specified domain
- A boundary face mask is applied at the domain midpoint in order fix the pressure value there
- The solver uses damped Newton iteration with initial damping factor of 0.05
- Surface charges should be specified in SI units (C/m²)
- Domain should be specified in SI units (meters)
"""
function mpbpsolve(;
        n = 21,
        chargenumbers = [-1, 1],
        domain = [0, 1.0e-10],
        surfacecharge = [0.16, -0.16],
        bulkmolarity = 0.1
    )

    X, grid = prepare_grid(n, domain)
    data = ICMPBP.ICMPBData(; z = chargenumbers, q = surfacecharge)
    ICMPBP.set_molarity!(data, bulkmolarity)
    sys = ICMPBP.ICMPBSystem(grid, data)
    sol = solve(sys; inival = 0.01, damp_initial = 0.05, verbose = "n")
    c0 = ICMPBP.calc_c0mol(sol, sys)
    c = ICMPBP.calc_cmol(sol, sys)
    return X, c0, c[1, :], c[2, :]
end

"""
    icmpbpsolve(; n=21, 
                  chargenumbers=[-1, 1], 
                  domain=[0, 1.0e-10], 
                  surfacecharge=[0.16, -0.16], 
                  averagemolarity=1)

Solve the ion-conserving modified Poisson-Boltzmann problem with specified average molarity.

This function solves the modified Poisson-Boltzmann equation with ion conservation constraints.
Unlike `mpbpsolve`, this solver conserves the total number of ions in the system and uses
an average molarity parameter rather than bulk molarity. This is particularly useful for
systems where ion conservation is physically important.

# Arguments
- `n::Int=21`: Number of grid points along the domain. Must be odd (automatically incremented if even). Default is 21.
- `chargenumbers::Vector{Int}=[-1, 1]`: Charge numbers of the ionic species (e.g., [-1, 1] for monovalent anions and cations).
- `domain::Vector{Float64}=[0, 1.0e-10]`: Domain boundaries in meters [start, end]. Default is [0, 1.0e-10] (0 to 0.1 nm).
- `surfacecharge::Vector{Float64}=[0.16, -0.16]`: Surface charge densities at the domain boundaries in C/m².
- `averagemolarity::Float64=1`: Average molarity of the electrolyte solution in mol/L.

# Returns
- `X`: Grid coordinates in m (Vector{Float64})
- `c0`: Solvent concentration in mol/L  (Vector{Float64})
- `c1`: Concentration of first ionic species along the grid in mol/L (Vector{Float64})
- `c2`: Concentration of second ionic species along the grid in mol/L (Vector{Float64})

# Examples
```julia
# Solve with default parameters
X, c0, c_anion, c_cation = icmpbpsolve()

# Solve with custom parameters and higher resolution
X, c0, c_anion, c_cation = icmpbpsolve(
    n=101,
    chargenumbers=[-2, 1],  # divalent anions, monovalent cations
    domain=[0, 5.0e-10],
    surfacecharge=[0.3, -0.3],
    averagemolarity=0.5
)
```

# Notes
- The function enforces an odd number of grid points for numerical stability
- Ion conservation is enforced through the `conserveions=true` parameter in ICMPBData
- The solver employs damped Newton iteration with initial damping factor of 0.05
- Average molarity represents the spatial average concentration, which may differ from bulk concentration due to ion conservation
- Surface charges should be specified in SI units (C/m²)
- Domain should be specified in SI units (meters)

# Differences from mpbpsolve
- Enforces ion conservation constraints
- Uses average rather than bulk molarity
- Requires odd number of grid points
- May show different concentration profiles due to conservation effects
"""
function icmpbpsolve(;
        n = 21,
        chargenumbers = [-1, 1],
        domain = [0, 1.0e-10],
        surfacecharge = [0.16, -0.16],
        averagemolarity = 1
    )

    X, grid = prepare_grid(n, domain)
    data = ICMPBP.ICMPBData(; z = chargenumbers, q = surfacecharge, conserveions = true)
    ICMPBP.set_molarity!(data, averagemolarity)
    sys = ICMPBP.ICMPBSystem(grid, data)
    sol = solve(sys; inival = unknowns(sys, data), damp_initial = 0.05, verbose = "n")
    c0 = ICMPBP.calc_c0mol(sol, sys)
    c = ICMPBP.calc_cmol(sol, sys)
    return X, c0, c[1, :], c[2, :]
end

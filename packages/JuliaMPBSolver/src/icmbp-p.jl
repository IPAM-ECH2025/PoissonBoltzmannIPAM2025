### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
if isdefined(Main, :PlutoRunner)
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using VoronoiFVM
    using ExtendableGrids
    using LinearAlgebra
    using LessUnitful
end

# ╔═╡ ab7a4a08-1b71-40e4-b940-8ae3bc441bfd
md"""
This file is part of the package and uses the fact that Pluto notebooks are 
just Julia scripts. It can be edited in any editor (please don't touch the UUIDs
describing the cell boundaries) or in the browser as any other Pluto notebook.
"""

# ╔═╡ ef660f6f-9de3-4896-a65e-13c60df5de1e
md"""
## Ion conserving MPB solver with Pressure
"""

# ╔═╡ 920b7d84-56c6-4958-aed9-fc67ba0c43f6
md"""
## Intro

This code implements the model described in
[Müller, R., Fuhrmann, J., & Landstorfer, M. (2020). Modeling polycrystalline electrode-electrolyte interfaces: The differential capacitance. Journal of The Electrochemical Society, 167(10), 106512](https://iopscience.iop.org/article/10.1149/1945-7111/ab9cca/meta)
in a symmetric cell with charged surfaces and ion conservation.


Equation numbers refer to the paper.

Concentrations are given in number densities.
"""

# ╔═╡ 87ac16f4-a4fc-4205-8fb9-e5459517e1b8
md"""
If not stated otherwise, all calculations and calculation results are in coherent SI units.
"""

# ╔═╡ 7d77ad32-3df6-4243-8bad-b8df4126e6ea
md"""
## Model data
"""

# ╔═╡ 4cabef42-d9f9-43fe-988e-7b54462dc775
md"""
#### ICMPBData
"""

# ╔═╡ 55b2ee36-c4f9-4ba3-84ed-faeb556aa026
makeδ(v, χ, T) = sqrt(v * χ * 3 * ph"ε_0" * ph"k_B" * T)

# ╔═╡ 0d825f88-cd67-4368-90b3-29f316b72e6e
begin
    """
        ICMPBData

    Data structure containing data for equilibrium calculations.
    All data including molarity in SI basic units
    """
    Base.@kwdef mutable struct ICMPBData

        "Ion charge numbers."
        z::Vector{Int} = [-1, 1]

        "Number of ionic species"
        N::Int64 = length(z)

        "Ion solvation numbers"
        κ::Vector{Float64} = fill(10.0, N)

        "Bulk molarity transformed to number density"
        molarity::Float64 = 0.1 * ph"N_A" / ufac"dm^3"

        "Bulk ion number densities"
        n_E::Vector{Float64} = fill(molarity, N)

        "Average ion number densities"
        n_avg::Vector{Float64} = fill(molarity, N)

        "Surface charges"
        q::Vector{Float64} = [0, 0] * ufac"C/m^2"

        "Solvent molarity"
        n0_ref::Float64 = 55.508 * ph"N_A" / ufac"dm^3"

        "Solvent molecular volume"
        v0::Float64 = 1 / n0_ref

        "Unsolvated ion molecular volume"
        vu::Vector{Float64} = fill(1 / n0_ref, N)

        "Dielectric susceptibility"
        χ0::Float64 = 78.49 - 1

        "Dielectric susceptibility"
        χ = fill(0.0, N)


        "Solvent molar fraction index"
        i0::Int = N + 1

        "Electric potential species index"
        iφ::Int = i0 + 1

        "Pressure species index"
        ip::Int = iφ + 1

        "Field strength species index"
        iE::Int = ip + 1

        "Offset of n_E in species list"
        coffset::Int = iE

        "Reference pressure"
        p_ref::Float64 = 1.0e5 * ufac"Pa"

        "Pressure scaling nparameter"
        pscale::Float64 = 1.0 * ufac"GPa"

        "Concentration scaling parameter"
        cscale::Float64 = ph"N_A"

        "Reference voltage"
        E_ref::Float64 = 0.0 * ufac"V"

        "Temperature"
        T::Float64 = 298.15 * ufac"K"

        δ0 = makeδ(v0, χ0, T)

        δ = [makeδ(κ[i] * v0 + vu[i], χ[i], T) for i in 1:N]


        "Temperature times Boltzmann constant"
        kT::Float64 = ph"k_B" * T

        "Electron charge"
        e::Float64 = ph"e"

        "Vacuum permittivity"
        ε_0::Float64 = ph"ε_0"

        "Ion conservation flag"
        conserveions::Bool = false

        "node volumes" # we should be able to query this from the system
        nv::Vector{Float64} = Float64[]

    end
end

# ╔═╡ 006ebe22-7bed-45db-acc4-d3e46f5ed7b7
ICMPBData()

# ╔═╡ 858ed8e1-84b1-4105-8ea0-45209aea40c6
md"""
#### `apply_charge!(data,q)`
"""

# ╔═╡ 4929c105-4c01-4c83-ad2f-2056a8c51d29
function apply_charge!(data::ICMPBData, q)
    data.q .= [q, - q]
    return data
end

# ╔═╡ f3049938-2637-401d-9411-4d7be07c19ca
md"""
#### set_molarity!(data,M)
"""

# ╔═╡ e69e10cc-e21a-418d-90b2-ae218dca0c73
md"""
Set the molarity of the electrolyte and update depending data
"""

# ╔═╡ 5d6340c4-2ddd-429b-a60b-3de5570a7398
function set_molarity!(data::ICMPBData, M_E)
    n_E = M_E * ph"N_A" / ufac"dm^3"
    data.molarity = n_E
    data.n_E = fill(n_E, data.N)
    data.n_avg = fill(n_E, data.N)
    return data
end

# ╔═╡ 30c6a176-935b-423f-9447-86f78746322f
md"""
#### debyelength(data)

```math
L_{Debye}=\sqrt{ \frac{(1+χ)ε_0k_BT}{e^2n_E}}
```
"""

# ╔═╡ a41c6c1f-ceb5-4590-a421-cae5078d167b
function L_Debye(data)
    return sqrt(
        (1 + data.χ) * data.ε_0 * ph"k_B" * data.T / (data.e^2 * data.n_E[1]),
    )
end;

# ╔═╡ a21545da-3b53-47af-b0c4-f253b37dc84f
md"""

#### dlcap0(data)
Double layer capacitance at $φ=0$
```math
C_{dl,0}=\sqrt{\frac{2(1+χ) ε_0e^2 n_E}{k_BT}}
```
"""

# ╔═╡ 1d22b09e-99c1-4026-9505-07bdffc98582
function dlcap0(data::ICMPBData)
    return sqrt(
        2 * (1 + data.χ) * ph"ε_0" * ph"e"^2 * data.n_E[1] / (ph"k_B" * data.T),
    )
end;

# ╔═╡ fe704fb4-d07c-4591-b834-d6cf2f4f7075
# ╠═╡ skip_as_script = true
#=╠═╡
let
    data=ICMPBData()
    set_molarity!(data,0.01)
    data.χ=78.49-1
    cdl0=dlcap0(data)
    @assert cdl0 ≈ 22.84669184882525ufac"μF/cm^2"
end
  ╠═╡ =#

# ╔═╡ 5a210961-19fc-40be-a5f6-033a80f1414d
md"""
Check with Bard/Faulkner: the value must be $(22.8)μF/cm^2")
"""

# ╔═╡ 5eca37ba-f858-45fb-a66a-3795327dfd18
md"""
## Model equations
"""

# ╔═╡ 84ed772f-f04a-4746-a864-3c9c6cc8bb67
function W(x)
    # W(x) = (coth(x) - 1/x) / x
    if abs(x) < 0.5
        u = x * x
        # Horner scheme with 4 terms thx 深度求索
        result = -1.0 / 4725.0    # u^3 coefficient
        result = 2.0 / 945.0 + u * result
        result = -1.0 / 45.0 + u * result
        result = 1.0 / 3.0 + u * result
        return 3 * result
    else
        return 3 * ((coth(x) - 1.0 / x) / x)
    end
end

# ╔═╡ 7ca19f45-b95e-4e21-a5aa-d4140669b3b8
function Λ(x) # thx 深度求索
    if abs(x) < 0.1
        u = x * x
        return log(1.0 + u * (1.0 / 6.0 + u * (1.0 / 120.0 + u * (1.0 / 5040.0 + u * (1.0 / 362880.0)))))
    else
        return log(sinh(x) / x)
    end
end

# ╔═╡ a26cf11b-0ce1-4c1d-a64d-1917178ff676
md"""
### Mole fractions
Equilibrium expression for mole fractions (``α≥0``) (16)
```math
y_α(φ,p)=y_α^E\exp\left(\frac{-z_αe}{k_BT}(φ- φ^E)-\frac{v_α}{k_BT}(p-p^E)\right)
```
"""

# ╔═╡ cdd1d359-08fa-45a1-a857-e19f2adefcab
md"""
#### y_α(φ,p,α,data)

Ion molar fractions
"""

# ╔═╡ 188f67d8-2ae8-474c-8e58-68b8b4fde02e
function y_α(φ, p, α, data, ddata)
    η_φ = data.z[α] * data.e * (φ - data.E_ref)
    η_p = ddata.v[α] * (p * data.pscale - data.p_ref)
    return ddata.y_E[α] * exp(-(η_φ + η_p) / (data.kT))
end;

# ╔═╡ f70eed13-a6c2-4d54-9f30-113367afaf7d
md"""
#### y0(p,data)

Solvent molar fraction
"""

# ╔═╡ d7531d5f-fc2d-42b2-9cf9-6a737b0f0f8d
function y0(p, data, ddata)
    return ddata.y0_E * exp(-data.v0 * (p * data.pscale - data.p_ref) / (data.kT))
end;

# ╔═╡ f6f004a6-d71b-4813-a363-9f51dc37e42a
md"""
### Poisson equation
Poisson equation (32a)

```math
-∇⋅(1+χ)ε_0∇φ = q(φ,p)
```
"""

# ╔═╡ 824c610b-6e5e-48a3-be37-19104f52d1d9
md"""
#### Space charge expression
"""

# ╔═╡ 2e11ce81-7d0d-498f-9ddd-7d4d836ab42f
md"""
Solvated ion volumes:
```math
	v_α=(1+κ_α)v_0
```
"""

# ╔═╡ b1e062c6-f245-4edc-aa02-871e2c776998
md"""
Incompressibility condition (14)
```math
\begin{aligned}
1&=∑\limits_α v_αn_α= n ∑\limits_α v_α y_α\\
n&=\frac1{∑\limits_α v_α y_α}
\end{aligned}
```
"""

# ╔═╡ c4cc940c-74aa-45f8-a2fa-6016d7c3c145
md"""
Space charge
```math
\begin{aligned}
q(φ,p)&=e∑\limits_α z_αn_α = ne∑\limits_α z_αy_α\\
      &=e\frac{∑\limits_α z_αy_α(\phi,p)}{∑\limits_α v_α y_α(\phi,p)}
\end{aligned}
```
"""

# ╔═╡ b07246b8-aec5-4161-8879-8cefb350aced
function spacecharge(u, data)
    (; iφ, ip, i0) = data
    y = u[i0]
    sumyz = zero(eltype(u))
    sumyv = data.v0 * y
    for α in 1:(data.N)
        y = u[α]
        sumyz += data.z[α] * y
        v = data.vu[α] + data.κ[α] * data.v0
        sumyv += v * y
    end
    return data.e * sumyz / sumyv
end

# ╔═╡ bf7f8bab-4807-440b-8796-fcf75ad313d7
function susceptibility(u, data)
    (; iE, i0) = data
    y = u[i0]
    χ = data.v0 * u[i0] * data.χ0*W(data.δ0 * u[iE] / data.kT)
    sumyv = data.v0 * y
    for α in 1:(data.N)
        v = data.vu[α] + data.κ[α] * data.v0
        χ += v * u[α] * data.χ[α]*W(data.δ[α] * u[iE] / data.kT)
        sumyv += v * u[α]
    end
    return χ / sumyv
end

# ╔═╡ b41838bb-3d5b-499c-9eb5-137c252ae366
md"""
#### Sum of mole fractions
"""

# ╔═╡ a468f43a-aa20-45dc-9c21-77f5adf2d700
function local_ysum(u, data)
    (; iφ, ip, i0) = data

    sumy = u[i0]
    for α in 1:(data.N)
        sumy += u[α]
    end
    return sumy
end

# ╔═╡ 13fc2859-496e-4f6e-8b22-36d9d55768b8
md"""
#### Derived data

Update derived data in data record.

Calculate bulk mole fractions from incompressibiltiy:
```math
\begin{aligned}
∑\limits_αv_αn_α^E&=1\\
n_0^E&=\frac1{v_0}\left(1-∑\limits_{α>0}v_αn_α^E\right)\\
n^E&=\frac1{v_0}\left(1-∑\limits_{α>0}v_αn_α^E\right)+ ∑\limits_{α>0}n_α^E\\
   &=\frac1{v_0}\left(1-∑\limits_{α>0}(v_α-v_0)n_α^E\right)\\
   &=\frac1{v_0}\left(1-∑\limits_{α>0}((1+ κ_α)v_0-v_0)n_α^E\right)\\
   &=\frac1{v_0}\left(1-∑\limits_{α>0}κ_αv_0n_α^E\right)\\
   &=\frac1{v_0}-∑\limits_{α>0}κ_αn_α^E\\
y_α^E&=\frac{n_α^E}{n^E}
\end{aligned}
```
"""

# ╔═╡ 7c5e9686-9d66-4ff0-84a7-c7b22596ab57
begin
    struct DerivedData{T}
        "Effective ion volumes"
        v::Vector{Float64}
        "Bulk ion mole fractions"
        y_E::Vector{T}
        "Bulk solvent mole fraction"
        y0_E::T
    end

    function DerivedData(data::ICMPBData, n_E)
        (; κ, v0, vu, T) = data
        c0 = zero(eltype(n_E)) + 1 / v0
        barc = zero(eltype(n_E))
        v = vu + κ * v0
        N = length(κ)
        for α in 1:N
            barc += n_E[α]
            c0 -= n_E[α] * (1 + κ[α])
        end
        barc += c0
        y_E = n_E / barc
        y0_E = c0 / barc
        return DerivedData(v, y_E, y0_E)
    end

    function DerivedData(data::ICMPBData)
        return DerivedData(data, data.n_E)
    end

end

# ╔═╡ b1e333c0-cdaa-4242-b71d-b54ff71aef83
let
    data = ICMPBData()
    set_molarity!(data, 0.01)
    ddata = DerivedData(data)
    sumyz = 0.0
    sumyv = ddata.y0_E * data.v0
    sumy = ddata.y0_E
    for α in 1:data.N
        v = (1.0 + data.κ[α]) * data.v0
        sumyz += ddata.y_E[α] * data.z[α]
        sumyv += ddata.y_E[α] * v
        sumy += ddata.y_E[α]
    end
    @assert sumy ≈ 1.0
end

# ╔═╡ 55bd7b9a-a191-4a0b-9c6b-13733be5023e
md"""
#### c_num!(c,φ,p, data)
Calculate number concentrations at discretization node
```math
\begin{aligned}
	n&=\sum_{i=0}^N y_α v_α\\
	n_α&=ny_α
\end{aligned}
```
"""

# ╔═╡ 3ceda3b1-bf1c-4126-b94f-2ee03e8dde99
function c_num!(c, y, data, ddata)
    (; i0) = data
    sumyv = data.v0 * y[i0]
    for α in 1:(data.N)
        c[α] = y[α]
        sumyv += y[α] * ddata.v[α]
    end
    return c ./= sumyv
end;

# ╔═╡ 97c5942c-8eb4-4b5c-8951-87ac0c9f396d
function c0_num(y, data, ddata)
    (; i0) = data

    y0 = y[i0]
    sumyv = data.v0 * y0
    for α in 1:(data.N)
        sumyv += y[α] * ddata.v[α]
    end
    return y0 / sumyv
end;

# ╔═╡ 0c54efd0-f279-4dc6-8b00-ba092dd13f44
md"""
#### calc_cnum(sol,sys)

Obtain ion number densities from system
"""

# ╔═╡ 800dfed8-9f29-4138-96f8-e8bf1f2f00e6
function calc_cnum(sol, sys)
    data = sys.physics.data
    i3 = sys.grid[BFaceNodes][3][1]
    if data.conserveions
        ddata = DerivedData(data, sol[(data.coffset + 1):end, i3])
    else
        ddata = DerivedData(data)
    end
    (; iφ, ip, N) = data
    grid = sys.grid
    nnodes = num_nodes(grid)
    conc = zeros(data.N, nnodes)
    for i in 1:nnodes
        @views c_num!(conc[:, i], sol[1:(N + 1), i], data, ddata)
    end
    return conc
end;

# ╔═╡ 24910762-7d56-446b-a758-d8e830fe9a09
function calc_c0num(sol, sys)
    data = sys.physics.data
    grid = sys.grid
    i3 = sys.grid[BFaceNodes][3][1]
    if data.conserveions
        ddata = DerivedData(data, sol[(data.coffset + 1):end, i3])
    else
        ddata = DerivedData(data)
    end
    (; iφ, ip, N) = data
    nnodes = num_nodes(grid)
    c0 = zeros(nnodes)
    for i in 1:nnodes
        @views c0[i] = c0_num(sol[1:(N + 1), i], data, ddata)
    end
    return c0
end;

# ╔═╡ 3b7a90cd-8f58-4abc-805a-2891ad6ceb9a
function calc_χ(sol, sys)
    data = sys.physics.data
    grid = sys.grid
    nnodes = num_nodes(grid)
    χ = zeros(nnodes)

    for i in 1:nnodes
        @views χ[i] = susceptibility(sol[:, i], data)
    end
    return χ
end

# ╔═╡ 9fe3ca93-c051-426e-8b9a-cc59f59319ad
md"""
#### calc_cmol(sol,sys)

Obtain ion  molarities (molar densities in mol/L)  from system
"""

# ╔═╡ 2ee34d76-7238-46c2-94d1-a40d8b017af6
calc_cmol(sol, sys) = calc_cnum(sol, sys) / (ph"N_A" * ufac"mol/dm^3");

# ╔═╡ 79cc671b-ef6e-42da-8641-61e43f221cb1
calc_c0mol(sol, sys) = calc_c0num(sol, sys) / (ph"N_A" * ufac"mol/dm^3");

# ╔═╡ 02d69f1c-4525-4f69-9938-cb0495171c3a
md"""
#### ysum(sys,sol)
"""

# ╔═╡ 48670f54-d303-4c3a-a191-06e6592a2e0a
function ysum(sys, sol)
    data = sys.physics.data
    n = size(sol, 2)
    sumy = zeros(n)
    for i in 1:n
        sumy[i] = local_ysum(sol[:, i], data)
    end
    return sumy
end

# ╔═╡ 7a607454-7b75-4313-920a-2dbdad258015
md"""
## The pressure equation
"""

# ╔═╡ 9cb8324c-896f-40f8-baa8-b7d47a93e9f5
md"""
This possibility to handle the pressure has been introduced in 

[J. Fuhrmann, “Comparison and numerical treatment of generalised Nernst–Planck models,” Computer Physics Communications, vol. 196, pp. 166–178, 2015.](https://dx.doi.org/10.1016/j.cpc.2015.06.004).

Starting with the momentum balance in mechanical equilibrium
```math
	\nabla p = -q\nabla \varphi
```
by taking the divergence on both sides of the equation, one derives the pressure Poisson problem
```math
\begin{aligned}
	-\Delta p &= \nabla\cdot q\nabla \varphi & \text{in}\; \Omega\\
      p&=p_{bulk} & \text{on}\; \Gamma_{bulk}\\
	(\nabla p + q\nabla \varphi)\cdot \vec n &=0 & \text{on}\; \partial\Omega\setminus\Gamma_{bulk}\\
\end{aligned}
```
"""

# ╔═╡ 003a5c0b-17c7-4407-ad23-21c0ac000fd4
md"""
The bulk Dirichlet boundary condition for the pressure is necessary to make the solution unique. It is reasonable to set the ``\varphi`` to a bulk value at ``\Gamma_{bulk}`` as well, and to calculate ``p_{bulk}`` from the molar fraction sum constraint.
"""

# ╔═╡ dc04ffc5-c96d-48c2-b89a-9094c57f1623
md"""
## Implementation
"""

# ╔═╡ 05695820-fa21-49b7-b52f-8a94cf2fa0fa
md"""
We use N+2 fields of unknowns in the following sequence:
``y_1 \dots y_N, y_0, \varphi, p``. Pressures are scaled by `pscale` (default: 1GPa).
"""

# ╔═╡ b9b0cb4f-cf72-418e-a65e-0f4c8a10e34c
md"""
#### `reaction!(f,u,node, data)`

Callback which runs in every grid point.
- Calculate space charge density and add this to the poisson equation
- Calculate molar fractions from potential and pressure if ion conservation is not required
"""

# ╔═╡ e1c13f1e-5b67-464b-967b-25e3a93e33d9
function reaction!(f, u, node, data)
    (; i0, iφ, ip, N) = data
    φ = u[iφ]
    p = u[ip]
    f[iφ] = -spacecharge(u, data)
    return
end;

# ╔═╡ c1168002-a716-4568-9a52-ac00f32f0134
md"""
#### `poisson_and_p_flux!(f, u, edge, data)`

Runs on every grid edge. Calculate fluxes for the Poisson and the pressure equations.
"""

# ╔═╡ 64e47917-9c61-4d64-a6a1-c6e8c7b28c59
function poisson_and_p_flux!(f, u, edge, data)
    (; iφ, ip, N) = data
    f[iφ] = (1.0 + data.χ0) * data.ε_0 * (u[iφ, 1] - u[iφ, 2])
    uu = zeros(eltype(u), N + 3)
    for i in 1:(N + 3)
        uu[i] = u[i, 1]
    end
    q1 = spacecharge(uu, data)
    for i in 1:(N + 3)
        uu[i] = u[i, 2]
    end
    q2 = spacecharge(uu, data)
    f[ip] =
        (u[ip, 1] - u[ip, 2]) + (u[iφ, 1] - u[iφ, 2]) * (q1 + q2) / (2 * data.pscale)
    return
end;

# ╔═╡ 05acd04f-74af-42f9-b039-7ee5b2ba63ff
md"""
#### `bcondition!(y, u, bnode, data)`

Boundary condition callback. The Dirichlet condition for the pressure in the mid of the domain ensures uniqueness of the pressure equation. 
"""

# ╔═╡ 743b9a7a-d6ac-4da0-8538-2045d965b547
function bcondition!(y, u, bnode, data)
    (; iφ, ip) = data
    boundary_neumann!(y, u, bnode, species = iφ, region = 2, value = data.q[2])
    boundary_neumann!(y, u, bnode, species = iφ, region = 1, value = data.q[1])
    boundary_dirichlet!(y, u, bnode, species = ip, region = 3, value = 0)
    return nothing
end

# ╔═╡ 1d4cd5e2-e0b7-4366-b7d0-dcbe65be738e
md"""
#### ionconservation!(f, u, sys, data)
"Generic callback" which shall ensure the ion conservation constraint.
This method runs over the full grid, and its sparsity pattern is automatically detected. It is called if ion conservation is required. In this case, ``n^E_\alpha`` are additional unknowns scaled by `cscale` (default: ``N_A``) which are attached to the mid of the domain (node `i3`), and additional equations need to be assembled. These are `N-1` ion conservation constraints and the an electroneutrality constraint.
"""

# ╔═╡ 2b779474-fe3b-4367-8553-3851341d719b
md"""
#### ICMPBSystem(grid, data)

Create MPB system with pressure. Ion conserving if `data.conserveions==true`. In that case, the solution is a sparse matrix.
"""

# ╔═╡ 50a6c3a0-2d24-4614-99e8-3645c8e7d5ba
md"""
#### unknowns(sys, data)
Initialize and return unknown vector.
"""

# ╔═╡ b0a45e53-8b98-4e18-8b41-7f6d0bc1f76e
function VoronoiFVM.unknowns(sys, data::ICMPBData)
    (; i0, iφ, ip, iE, coffset, N) = data
    u = unknowns(sys, inival = 0)
    i3 = sys.grid[BFaceNodes][3][1]
    for α in 1:N
        u[α, :] .= 0.1
        if data.conserveions
            u[coffset + α, i3] = data.n_E[α] / data.cscale
        end
    end
    u[i0, :] .= 1 - N * 0.1
    return u
end

# ╔═╡ dbccaa88-65d9-47ab-be78-83df64a6db24
function ionconservation!(f, u, sys, data)
    (; coffset, i0, iφ, iE, ip, N, z, nv, n_avg) = data
    # Set the result to zero
    f .= 0
    # Find  mid-of-the-domain node number from boundary region 3
    i3 = sys.grid[BFaceNodes][3][1]
    X = sys.grid[Coordinates][1, :]

    # Parameters u and f come as vectors, `idx` allows to access their contents with
    # two-dimensional indexing. This might be changed in a later Version of VoronoiFVM.
    idx = unknown_indices(unknowns(sys))

    # Obtain values of the bulk molecular densities
    if data.conserveions
        n_E = [u[idx[coffset + i, i3]] * data.cscale for i in 1:N]
        # Calculate derived data
        ddata = DerivedData(data, n_E)
    else
        ddata = DerivedData(data)
    end

    # Calculate molar fractions for all nodes of the grid
    for i in 1:num_nodes(sys.grid)
        f[idx[i0, i]] = u[idx[i0, i]] - y0(u[idx[ip, i]], data, ddata)
        for α in 1:N
            f[idx[α, i]] = u[idx[α, i]] - y_α(u[idx[iφ, i]], u[idx[ip, i]], α, data, ddata)
        end
        if i == 1
            f[idx[iE, i]] = u[idx[iE, i]] - abs((u[idx[iφ, i + 1]] - u[idx[iφ, i]]) / (X[i + 1] - X[i]))
        elseif i == num_nodes(sys.grid)
            f[idx[iE, i]] = u[idx[iE, i]] - abs((u[idx[iφ, i - 1]] - u[idx[iφ, i]]) / (X[i - 1] - X[i]))
        else
            f[idx[iE, i]] = u[idx[iE, i]] - abs((u[idx[iφ, i + 1]] - u[idx[iφ, i - 1]]) / (X[i + 1] - X[i - 1]))
        end
    end

    if data.conserveions
        # Get size of the domain
        L = sum(nv)
        # Initialize electroneutrality constraint for n^E_N
        f[idx[coffset + N, i3]] = u[idx[coffset + N, i3]]
        for α in 1:(N - 1)
            # Initialize ion conservation constrain for n_α
            f[idx[coffset + α, i3]] = -n_avg[α] * L / data.cscale

            # Update electroneutrality constraint for n^E_N
            f[idx[coffset + N, i3]] += z[α] * u[idx[coffset + α, i3]] / z[N]
        end

        y = zeros(eltype(u), N + 1)
        uu = zeros(eltype(u), N + 1)
        # Calculate number density integrals
        for iv in 1:length(nv)
            for α in 1:(N + 1)
                uu[α] = u[idx[α, iv]]
            end
            c_num!(y, uu, data, ddata)

            for α in 1:(N - 1)
                # Update ion conservation constraint for n_α
                f[idx[coffset + α, i3]] += y[α] * nv[iv] / data.cscale
            end
        end
    end
    return nothing
end

# ╔═╡ 7bf3a130-3b47-428e-916f-4a0ec1237844
function ICMPBSystem(grid, data)

    data.nv = ones(num_nodes(grid)) # trigger sparsity detector
    sys = VoronoiFVM.System(
        grid;
        data = data,
        flux = poisson_and_p_flux!,
        reaction = reaction!,
        bcondition = bcondition!,
        generic = ionconservation!,
        unknown_storage = :sparse
    )

    # Enable species for all fields
    for i in 1:data.N
        enable_species!(sys, i, [1])
    end

    enable_species!(sys, data.i0, [1])
    enable_species!(sys, data.iφ, [1])
    enable_species!(sys, data.ip, [1])
    enable_species!(sys, data.iE, [1])

    # Enable species in mid of the domain for ion conservation and
    # electroneutrality constraint
    if data.conserveions
        for α in 1:data.N
            enable_boundary_species!(sys, data.coffset + α, [3])
        end
        data.nv = nodevolumes(sys)
    end

    return sys
end;

# ╔═╡ 3dc0d408-ab00-4da0-9999-2ddf6e4fbf60
md"""
## Capacitance calculation
"""

# ╔═╡ a67c0d46-456d-4b0c-8519-961b37043350
md"""
#### qsweep(sys)

Sweep over series of surface charges and calculate resulting potential
difference.
"""

# ╔═╡ 3b389ecf-4c63-4eb8-b8be-76442eacef80


# ╔═╡ 178b947f-3fef-44ed-9eca-fdb9916bc2b6
function qsweep(sys; qmax = 10, nsteps = 100, verbose = "", kwargs...)
    data = deepcopy(sys.physics.data)
    (; ip, iφ) = data
    apply_charge!(data, 0 * ph"e" / ufac"nm^2")
    state = VoronoiFVM.SystemState(sys; data)
    sol = solve!(state; inival = unknowns(sys, data), damp_initial = 0.1, verbose, kwargs...)

    volts = []
    Q = []
    for q in range(0, qmax, length = 50)
        apply_charge!(data, q * ph"e" / ufac"nm^2")
        sol = solve!(state; inival = sol, damp_initial = 0.1, verbose, kwargs...)
        push!(volts, (sol[iφ, end] - sol[iφ, 1]) / 2)
        # Division by 2 comes in because the voltage we get here is the difference
        # between the electrodes and not the difference between electrode and bulk
        # which would correspond to the usual half-cell dlcap experiment
        push!(Q, q * ph"e" / ufac"nm^2")
    end
    dlcaps = -(Q[2:end] - Q[1:(end - 1)]) ./ (volts[2:end] - volts[1:(end - 1)])
    return volts[1:(end - 1)], dlcaps
end

# ╔═╡ 4032bc46-7820-45d4-bc11-9350ecf1797a
md"""
#### capscalc(sys, molarities)

Calculate double layer capacitances using qsweep results.

This provides an  "inverse" method to calculate these capacitances. Usually
one calculates charges dependent on voltages, here we calculate voltages dependent on charges.
"""

# ╔═╡ fc84996b-02c0-4c16-8632-79f4e1900f78
function capscalc(sys, molarities; kwargs...)
    result = []
    for imol in 1:length(molarities)
        data = sys.physics.data
        set_molarity!(data, molarities[imol])

        t = @elapsed volts, caps = qsweep(sys; kwargs...)
        cdl0 = dlcap0(data)
        push!(
            result,
            (
                voltages = volts,
                dlcaps = caps,
                cdl0 = cdl0,
                molarity = molarities[imol],
            ),
        )
    end
    return result
end

# ╔═╡ Cell order:
# ╟─ab7a4a08-1b71-40e4-b940-8ae3bc441bfd
# ╟─ef660f6f-9de3-4896-a65e-13c60df5de1e
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═920b7d84-56c6-4958-aed9-fc67ba0c43f6
# ╟─87ac16f4-a4fc-4205-8fb9-e5459517e1b8
# ╟─7d77ad32-3df6-4243-8bad-b8df4126e6ea
# ╟─4cabef42-d9f9-43fe-988e-7b54462dc775
# ╠═55b2ee36-c4f9-4ba3-84ed-faeb556aa026
# ╠═0d825f88-cd67-4368-90b3-29f316b72e6e
# ╠═006ebe22-7bed-45db-acc4-d3e46f5ed7b7
# ╟─858ed8e1-84b1-4105-8ea0-45209aea40c6
# ╠═4929c105-4c01-4c83-ad2f-2056a8c51d29
# ╟─f3049938-2637-401d-9411-4d7be07c19ca
# ╟─e69e10cc-e21a-418d-90b2-ae218dca0c73
# ╠═5d6340c4-2ddd-429b-a60b-3de5570a7398
# ╟─30c6a176-935b-423f-9447-86f78746322f
# ╠═a41c6c1f-ceb5-4590-a421-cae5078d167b
# ╟─a21545da-3b53-47af-b0c4-f253b37dc84f
# ╠═1d22b09e-99c1-4026-9505-07bdffc98582
# ╠═fe704fb4-d07c-4591-b834-d6cf2f4f7075
# ╟─5a210961-19fc-40be-a5f6-033a80f1414d
# ╟─5eca37ba-f858-45fb-a66a-3795327dfd18
# ╠═84ed772f-f04a-4746-a864-3c9c6cc8bb67
# ╠═7ca19f45-b95e-4e21-a5aa-d4140669b3b8
# ╟─a26cf11b-0ce1-4c1d-a64d-1917178ff676
# ╟─cdd1d359-08fa-45a1-a857-e19f2adefcab
# ╠═188f67d8-2ae8-474c-8e58-68b8b4fde02e
# ╟─f70eed13-a6c2-4d54-9f30-113367afaf7d
# ╠═d7531d5f-fc2d-42b2-9cf9-6a737b0f0f8d
# ╟─f6f004a6-d71b-4813-a363-9f51dc37e42a
# ╟─824c610b-6e5e-48a3-be37-19104f52d1d9
# ╟─2e11ce81-7d0d-498f-9ddd-7d4d836ab42f
# ╟─b1e062c6-f245-4edc-aa02-871e2c776998
# ╟─c4cc940c-74aa-45f8-a2fa-6016d7c3c145
# ╠═b07246b8-aec5-4161-8879-8cefb350aced
# ╠═bf7f8bab-4807-440b-8796-fcf75ad313d7
# ╟─b41838bb-3d5b-499c-9eb5-137c252ae366
# ╠═a468f43a-aa20-45dc-9c21-77f5adf2d700
# ╟─13fc2859-496e-4f6e-8b22-36d9d55768b8
# ╠═7c5e9686-9d66-4ff0-84a7-c7b22596ab57
# ╠═b1e333c0-cdaa-4242-b71d-b54ff71aef83
# ╟─55bd7b9a-a191-4a0b-9c6b-13733be5023e
# ╠═3ceda3b1-bf1c-4126-b94f-2ee03e8dde99
# ╠═97c5942c-8eb4-4b5c-8951-87ac0c9f396d
# ╟─0c54efd0-f279-4dc6-8b00-ba092dd13f44
# ╠═800dfed8-9f29-4138-96f8-e8bf1f2f00e6
# ╠═24910762-7d56-446b-a758-d8e830fe9a09
# ╠═3b7a90cd-8f58-4abc-805a-2891ad6ceb9a
# ╟─9fe3ca93-c051-426e-8b9a-cc59f59319ad
# ╠═2ee34d76-7238-46c2-94d1-a40d8b017af6
# ╠═79cc671b-ef6e-42da-8641-61e43f221cb1
# ╟─02d69f1c-4525-4f69-9938-cb0495171c3a
# ╠═48670f54-d303-4c3a-a191-06e6592a2e0a
# ╟─7a607454-7b75-4313-920a-2dbdad258015
# ╟─9cb8324c-896f-40f8-baa8-b7d47a93e9f5
# ╟─003a5c0b-17c7-4407-ad23-21c0ac000fd4
# ╟─dc04ffc5-c96d-48c2-b89a-9094c57f1623
# ╟─05695820-fa21-49b7-b52f-8a94cf2fa0fa
# ╟─b9b0cb4f-cf72-418e-a65e-0f4c8a10e34c
# ╠═e1c13f1e-5b67-464b-967b-25e3a93e33d9
# ╟─c1168002-a716-4568-9a52-ac00f32f0134
# ╠═64e47917-9c61-4d64-a6a1-c6e8c7b28c59
# ╟─05acd04f-74af-42f9-b039-7ee5b2ba63ff
# ╠═743b9a7a-d6ac-4da0-8538-2045d965b547
# ╟─1d4cd5e2-e0b7-4366-b7d0-dcbe65be738e
# ╠═dbccaa88-65d9-47ab-be78-83df64a6db24
# ╟─2b779474-fe3b-4367-8553-3851341d719b
# ╠═7bf3a130-3b47-428e-916f-4a0ec1237844
# ╟─50a6c3a0-2d24-4614-99e8-3645c8e7d5ba
# ╠═b0a45e53-8b98-4e18-8b41-7f6d0bc1f76e
# ╟─3dc0d408-ab00-4da0-9999-2ddf6e4fbf60
# ╟─a67c0d46-456d-4b0c-8519-961b37043350
# ╠═3b389ecf-4c63-4eb8-b8be-76442eacef80
# ╠═178b947f-3fef-44ed-9eca-fdb9916bc2b6
# ╟─4032bc46-7820-45d4-bc11-9350ecf1797a
# ╠═fc84996b-02c0-4c16-8632-79f4e1900f78

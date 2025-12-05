### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using Revise
    using PlutoUI, HypertextLiteral, UUIDs
    using VoronoiFVM
    using ExtendableGrids
    using LinearAlgebra
    using LessUnitful
    using Test
    using PythonPlot
    using JuliaMPBSolver.ICMPBP: ICMPBData, ICMPBSystem, L_Debye, set_molarity!, dlcap0, DerivedData, apply_charge!, ysum, qsweep, capscalc, calc_cmol, calc_c0mol, calc_χ, W
end

# ╔═╡ ef660f6f-9de3-4896-a65e-13c60df5de1e
md"""
# Ion conserving MPB with pressure, solvation, dielectric decrement
"""

# ╔═╡ 920b7d84-56c6-4958-aed9-fc67ba0c43f6
md"""
This code implements the model described in
[Müller, R., Fuhrmann, J., & Landstorfer, M. (2020). Modeling polycrystalline electrode-electrolyte interfaces: The differential capacitance. Journal of The Electrochemical Society, 2020, 167(10), 106512](https://iopscience.iop.org/article/10.1149/1945-7111/ab9cca/meta) with surface charge boundary conditions in a symmetric cell and ion conservation constraint.

For the later, see e.g. H. Sugioka, "Ion-Conserving Modified Poisson–Boltzmann Theory Considering a Steric Effect in an Electrolyte", Journal of the Physical Society of Japan, 2016, Vol. 85, No. 12, DOI [10.7566/JPSJ.85.124006](https://doi.org/10.7566/JPSJ.85.124006)


"""

# ╔═╡ aa9cfb26-ad55-4f37-804a-b33b706e9427
md"""
## Poisson equation
Domain: ``\Omega=(0,L)``
```math
-∇⋅(1+χ)ε_0∇φ = q(φ,p)
```

Space charge
```math
\begin{aligned}
q(φ,p)&=e∑\limits_α z_αn_α = ne∑\limits_{α=1}^N z_αy_α
      =e\frac{∑\limits_α z_αy_α(\phi,p)}{∑\limits_{α=0}^N v_α y_α(\phi,p)}
\end{aligned}
```
where 
```math
y_α(φ,p)=y_α^E\exp\left(\frac{-z_αe}{k_BT}(φ- φ^E)-\frac{v_α}{k_BT}(p-p^E)\right)
```

"""


# ╔═╡ 4d571068-133f-4b1a-86f6-3225439652eb
md"""
## Parameters
"""

# ╔═╡ 1facab2a-6ebf-492f-978f-2ace26b3bc4f
md"""
Here, for ``\alpha=1\dots N``, ``z_α`` are the charge numbers, and ``v_α=v_α^u + κ_α v_0`` are the effective molecular volumes, with solvation numbers ``κ_α`` and unsolvated  molecular volumes ``v_α^u``.
"""

# ╔═╡ 275536e9-394d-4686-b37c-ef10a623b400
md"""
The bulk mole fractions ``y_α^E``  are calculated
 from  the incompressibiltiy constraint 
`` ∑\limits_{α=1}^N  v_αn_α^E=1`` and given bulk number densities ``n_α^E`` (``α=1…N``):
```math
\begin{aligned}
n_0^E&=\frac1{v_0}\left(1-∑\limits_{α=1}^Nv_αn_α^E\right)\\
n^E&=\frac1{v_0}-∑\limits_{α=1}^Nκ_αn_α^E\\
y_α^E&=\frac{n_α^E}{n^E}
\end{aligned}
```
"""


# ╔═╡ e89e4da3-a2e4-41a3-99d9-327a9b862d28
md"""
## Surface with given charge `q`
```math
\begin{aligned}
	(1+χ)\varepsilon_0\nabla \phi|_{z=0} &= q\\
	-(1+χ)\varepsilon_0\nabla \phi|_{z=L} &= -q
\end{aligned}
```

"""

# ╔═╡ 63d99173-c779-4e72-b96f-15b27aec55cf
md"""
## Pressure equation
Pressure is calculated according to [J. Fuhrmann, “Comparison and numerical treatment of generalised Nernst–Planck models,” Computer Physics Communications, vol. 196, pp. 166–178, 2015.](https://dx.doi.org/10.1016/j.cpc.2015.06.004).

Starting with the momentum balance in mechanical equilibrium
```math
	\nabla p = -q\nabla \varphi
```
by taking the divergence on both sides of the equation, one derives the pressure Poisson problem
```math
\begin{aligned}
	-\Delta p &= \nabla\cdot q\nabla \varphi & \text{in}\; \Omega\\
	(\nabla p + q\nabla \varphi)\cdot \vec n &=0 & \text{on}\; \partial\Omega\\
\end{aligned}
```
In order to obtain uniqueness, the pressure in the center of the domain is set to zero:
```math
      p|_{x=\frac{L}2}=0
```
"""

# ╔═╡ bab7b779-339d-4702-a88a-e6cbf5a72cd0
md"""
## Ion conservation

Ion conservation is expressed by the constraint
```math
	\frac1{|\Omega|}\int\limits_\Omega n_α d\omega = n_α^{avg}\quad (α=1\dots N).
```
In order to implement this condition, the constraints (for ``α=1\dots N-1``) are added to the system, and the bulk molarities ``n_α^E`` are made variables with the electroneutrality constraint ``\sum\limits_{\alpha=1}^N z_α n_α^E =0``  .
"""

# ╔═╡ 16a79117-e89b-4478-a571-8c011b5784c1
md"""
## Simulation result
"""

# ╔═╡ 760e5861-7a6f-41bb-8aec-5e7466c6ec9f
md"""
## Simulation setup and run
"""

# ╔═╡ f4facb34-1f4a-432d-8a1e-30299e542bcd
begin
    const nm = ufac"nm"
    const V = ufac"V"
    const cm = ufac"cm"
    const dm = ufac"dm"
    const μF = ufac"μF"
    const N_A = ph"N_A"
    const mol = ufac"mol"
end;

# ╔═╡ ae11bded-9f67-4004-8786-ed54e1ccb932
surfcharge(n) = n * ph"e" / ufac"nm^2"

# ╔═╡ f75f1d3a-47e5-475b-97b1-bb275a510783
md"""
### Helper functions
"""

# ╔═╡ dc05f31c-a28e-4470-8916-72dda567b149
myround(x) = round(x, sigdigits = 4)

# ╔═╡ 000e87b5-cd1b-4b23-99be-6b7006502312
function capsplot(ax, result, title)
    hmol = 1 / length(result)
    ax.set_title(title)
    ax.grid()
    ax.set_xlabel("U/V")
    ax.set_ylabel(L"C_{dl}/(μF/cm^3)")

    for imol in 1:length(result)
        c = (1 - imol * hmol, 0, imol * hmol)

        ax.plot(
            result[imol].voltages, result[imol].dlcaps / (μF / cm^2),
            color = c, label = "$(result[imol].molarity)M"
        )
        ax.scatter(
            [0], [result[imol].cdl0] / (μF / cm^2),
            color = c
        )
    end
    return
end


# ╔═╡ 0a824568-997a-4a48-9fe3-71cc5f9b7471
begin
    function floataside(text::Markdown.MD; top = 1)
        uuid = uuid1()
        return @htl(
            """
            		<style>


            		@media (min-width: calc(700px + 30px + 300px)) {
            			aside.plutoui-aside-wrapper-$(uuid) {

            	color: var(--pluto-output-color);
            	position:fixed;
            	right: 1rem;
            	top: $(top)px;
            	width: 400px;
            	padding: 10px;
            	border: 3px solid rgba(0, 0, 0, 0.15);
            	border-radius: 10px;
            	box-shadow: 0 0 11px 0px #00000010;
            	/* That is, viewport minus top minus Live Docs */
            	max-height: calc(100vh - 5rem - 56px);
            	overflow: auto;
            	z-index: 40;
            	background-color: var(--main-bg-color);
            	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);

            			}
            			aside.plutoui-aside-wrapper > div {
            #				width: 300px;
            			}
            		}
            		</style>

            		<aside class="plutoui-aside-wrapper-$(uuid)">
            		<div>
            		$(text)
            		</div>
            		</aside>

            		"""
        )
    end
    floataside(stuff; kwargs...) = floataside(md"""$(stuff)"""; kwargs...)
end;


# ╔═╡ 2ef6b8d7-e5f3-4700-8a40-8feffab3569f
floataside(
    md"""
    - ``L/nm``: $(@bind L0 PlutoUI.Slider(2:1:10, default=10, show_value=true))		   
    - ``M_{avg}/(mol/dm^3)``:  $(@bind M1_avg PlutoUI.Slider(0.1:0.1:2, default=1, show_value=true))
    - ``n_e/(e/nm^2)``: $(@bind n1_e PlutoUI.Slider(0:0.5:5, default=1, show_value=true))
    - ``κ``: $(@bind kappa1 PlutoUI.Slider(0:1:20, default=10, show_value=true))
      ``M_{scale}/(mol/dm^3)``: $(@bind Mscale PlutoUI.Slider(5:5:60, default=60, show_value=true))
    """, top = 100
)


# ╔═╡ eacdd772-1869-406a-b601-64cdd6453ec1
begin
    data1 = ICMPBData(
        κ = [kappa1, kappa1]
    )
    set_molarity!(data1, M1_avg)
    data1.conserveions = true
    data1.χvar = true
end

# ╔═╡ f8c1c2bd-7466-491e-9132-4f15edcfa4c7
let
    clf()
    X = -10.0e8:1.0e6:10.0e8
    plot(X, W.(data1.δ0 * X / data1.kT))
    gcf()
end

# ╔═╡ a629e8a1-b1d7-42d8-8c17-43475785218e
begin
    L = L0 * nm
    n = 101
    X = range(0, L, length = n)
    grid = ExtendableGrids.simplexgrid(X)
    bfacemask!(grid, [L / 2], [L / 2], 3, tol = 1.0e-10 * nm)
    i3 = grid[BFaceNodes][3][1]
end

# ╔═╡ 8433319f-2f78-494c-9b2e-a5390cf93b00
sys1 = ICMPBSystem(grid, data1, valuetype = Float64);

# ╔═╡ 70910bd5-b8ca-4021-9b40-233b50ea5601
inival1 = unknowns(sys1, data1);

# ╔═╡ f579cf2d-9511-48a8-bf11-7400ef76ee3d
function plotsol(
        sol, sys; data = data(sys),
        grid = sys.grid, size = (600, 600), Mscale = 60
    )
    PythonPlot.clf()
    fig, ax = pyplot.subplots(3, 1)
    fig.set_size_inches(size[1] / 100, size[2] / 100)
    ax1 = ax[0]
    ax2 = ax[1]
    ax3 = ax[2]

    ax1.grid()
    i3 = sys.grid[BFaceNodes][3][1]
    ax1r = ax1.twinx()
    X = grid[Coordinates][1, :]
    #  ε_r = relpermittivity(sol, data; grid)
    c = calc_cmol(sol, sys)
    c0 = calc_c0mol(sol, sys)
    ax1.set_title("ϕ∈$(round.(Float64.(extrema(sol[data.iφ, :])), sigdigits = 3))")
    #, ε_r ∈$(round.(Float64.(extrema(ε_r)), sigdigits = 3))")
    ax1.plot(X / nm, sol[data.iφ, :], color = "green", linewidth = 2, label = "ϕ/V")
    ax1r.plot(X / nm, sol[data.ip, :], color = "pink", linewidth = 3, label = L"p/GPa")
    #   ax1.set_ylim(-10, 10)
    ax1.set_xlabel("z/nm")
    ax1.set_ylabel("ϕ/V")
    ax1.legend(loc = (0.1, 0.1))
    ax1r.legend(loc = (0.8, 0.1))
    #  ax1r.set_ylim(0, 80)

    ax2.grid()

    nv = nodevolumes(sys)
    cm, cp = c[1, :] ⋅ nv / L, c[2, :] ⋅ nv / L
    if data1.conserveions
        M_bulk = sol[7:8, i3] * data.cscale / (ph"N_A" * ufac"mol/dm^3")
        crm, crp = M_bulk[1], M_bulk[2]
        #  ax2.set_title("M_bulk=$(myround.((crm, crp))),  M_avg=$(myround.((cm, cp)))")
    end
    #  ax2.set_title("M_avg=$(myround.((cm, cp)))")
    ax2.set_xlabel("z/nm")
    ax2.set_ylabel("c/(mol/L)")
    ax2.set_ylim(0, Mscale)

    ax2.plot(X / nm, c[1, :], color = "blue", linewidth = 2, label = L"c^-")
    ax2.plot(X / nm, c[2, :], color = "red", linewidth = 2, label = L"c^+")
    ax2.plot(X / nm, c0, color = "green", linewidth = 2, label = L"c_{solvent}")
    ax2.legend(loc = (0.4, 0.1))


    ax3.grid()
    ax3.plot(X / nm, sol[data.iE, :], color = "pink", label = "E")
    ax3r = ax3.twinx()
    ax3r.plot(X / nm, calc_χ(sol, sys), color = "orange", label = "χ")
    ax3.legend(loc = (0.1, 0.1))
    ax3r.legend(loc = (0.8, 0.1))
    ax3r.set_ylim(0, 100)
    tight_layout()
    return PythonPlot.gcf()
end


# ╔═╡ d2df6ed0-e6f5-4677-b790-bfc40de7fd6a
begin
    Q = surfcharge(n1_e)
    sol1 = inival1
    for q in range(0, Q, length = 11)
        data1.q .= [-q, q]
        global sol1 = solve(sys1; inival = sol1, verbose = "")
    end
end

# ╔═╡ 5f153fe4-4476-401e-8674-b6902213c19e
md"""
Newton steps: $(length(history(sol1)))
"""

# ╔═╡ c7a08779-53f3-4fab-8bd4-3dffe6135c3b
plotsol(sol1, sys1; Mscale)

# ╔═╡ 0e4ec7f0-0aa8-4a32-96a3-40f63f32a12d
sol1

# ╔═╡ 7d7ebb45-2fb3-40ea-83f2-62d0a240b2db
@test isa(sol1, AbstractMatrix)

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─ef660f6f-9de3-4896-a65e-13c60df5de1e
# ╟─920b7d84-56c6-4958-aed9-fc67ba0c43f6
# ╟─aa9cfb26-ad55-4f37-804a-b33b706e9427
# ╟─4d571068-133f-4b1a-86f6-3225439652eb
# ╟─1facab2a-6ebf-492f-978f-2ace26b3bc4f
# ╟─275536e9-394d-4686-b37c-ef10a623b400
# ╟─e89e4da3-a2e4-41a3-99d9-327a9b862d28
# ╟─63d99173-c779-4e72-b96f-15b27aec55cf
# ╟─bab7b779-339d-4702-a88a-e6cbf5a72cd0
# ╟─16a79117-e89b-4478-a571-8c011b5784c1
# ╟─5f153fe4-4476-401e-8674-b6902213c19e
# ╠═c7a08779-53f3-4fab-8bd4-3dffe6135c3b
# ╠═0e4ec7f0-0aa8-4a32-96a3-40f63f32a12d
# ╠═eacdd772-1869-406a-b601-64cdd6453ec1
# ╟─760e5861-7a6f-41bb-8aec-5e7466c6ec9f
# ╠═f4facb34-1f4a-432d-8a1e-30299e542bcd
# ╟─2ef6b8d7-e5f3-4700-8a40-8feffab3569f
# ╠═a629e8a1-b1d7-42d8-8c17-43475785218e
# ╠═ae11bded-9f67-4004-8786-ed54e1ccb932
# ╠═8433319f-2f78-494c-9b2e-a5390cf93b00
# ╠═70910bd5-b8ca-4021-9b40-233b50ea5601
# ╠═d2df6ed0-e6f5-4677-b790-bfc40de7fd6a
# ╠═7d7ebb45-2fb3-40ea-83f2-62d0a240b2db
# ╠═f8c1c2bd-7466-491e-9132-4f15edcfa4c7
# ╟─f75f1d3a-47e5-475b-97b1-bb275a510783
# ╠═dc05f31c-a28e-4470-8916-72dda567b149
# ╠═f579cf2d-9511-48a8-bf11-7400ef76ee3d
# ╠═000e87b5-cd1b-4b23-99be-6b7006502312
# ╟─0a824568-997a-4a48-9fe3-71cc5f9b7471

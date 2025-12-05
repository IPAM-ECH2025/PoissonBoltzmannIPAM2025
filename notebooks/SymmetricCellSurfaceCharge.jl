### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ a70cef7d-2a2f-4155-bdf3-fec9df94c63f
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
    using DoubleFloats
    using PythonPlot
    using PythonCall
    using PythonPlot: pyplot
    using LaTeXStrings
    using Colors
    using DrWatson, PoissonBoltzmannIPAM2025
    using JuliaMPBSolver.ICMPBP: ICMPBP, ICMPBData, SurfaceChargedSymmetricCell, AbstractSymmetricCell, set_molarity!, calc_cmol, calc_c0mol, calc_χ, get_E, get_φ, get_p, get_c0,
        set_κ!, set_q!
end

# ╔═╡ af6ae00d-f032-4743-878b-e575466b6e84
md"""
# Symmetric cell with ion conservation and applied charge
"""

# ╔═╡ 87d66a6f-b40d-4d76-afc6-4b4f086e80a4
begin
    const nm = ufac"nm"
    const V = ufac"V"
    const cm = ufac"cm"
    const dm = ufac"dm"
    const μF = ufac"μF"
    const N_A = ph"N_A"
    const mol = ufac"mol"
end

# ╔═╡ 087614cc-a4f6-4867-a6bd-201d1f2a7fc2
begin
    const Ll = 10nm
    const Ls = 2nm
    const n = 101
end

# ╔═╡ c7a51cb9-790a-4818-95d7-c32f2252e8b3
Xl = range(0, Ll, length = n)

# ╔═╡ d87781e8-af89-41d0-a56d-94744406042d
Xs = range(0, Ls, length = n)

# ╔═╡ 4ef3d13c-2f57-4575-bc88-c8092ae6910f
begin
    gridl = simplexgrid(Xl)
    bfacemask!(gridl, [Ll / 2], [Ll / 2], 3, tol = 1.0e-10 * nm)

end

# ╔═╡ 4d78d81c-f1a3-4386-adce-2a7fa1e3eaef
begin
    grids = simplexgrid(Xs)
    bfacemask!(grids, [Ls / 2], [Ls / 2], 3, tol = 1.0e-10 * nm)

end

# ╔═╡ 3a99a9bb-7862-41f8-a535-d3c580da6909
md"""
All values are given with respect to SI basic units (m, kg, s, V, A)
"""

# ╔═╡ b24b7e23-61ea-41fc-a345-286e904c042b
data = ICMPBData(χvar = true, conserveions = true)

# ╔═╡ 1bb47749-edde-4bee-be9f-059a7652b354
begin
    smallcell = SurfaceChargedSymmetricCell(grids, data, dielectric_decrement = false)
    largecell = SurfaceChargedSymmetricCell(gridl, data, dielectric_decrement = false)

    smallcelldd = SurfaceChargedSymmetricCell(grids, data, dielectric_decrement = true)
    largecelldd = SurfaceChargedSymmetricCell(gridl, data, dielectric_decrement = true)

end;

# ╔═╡ 6aae7cca-01f1-45be-9d4b-a3942bd042f3
p = plotcells(largecell, largecelldd; figname = "largecell");p

# ╔═╡ 27c9e166-df00-4a30-a6ba-32e9a071136b
@test isa(p, Figure)

# ╔═╡ 773e6c0a-08ea-47b2-8924-42759f944c6b
p2 = plotcells(smallcell, smallcelldd, figname = "smallcell");p2

# ╔═╡ ceaff98b-1ba6-43b3-bd2f-0e7fe585f428
@test isa(p2, Figure)

# ╔═╡ 8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
html"""<hr>"""

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ afe4745f-f9f1-4e23-8735-cbec6fb79c41
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


# ╔═╡ b8fd36a7-d8d1-45f7-b66e-df9132168bfc
# https://discourse.julialang.org/t/adding-a-restart-process-button-in-pluto/76812/5
restart_button() = html"""
<script>
	const button = document.createElement("button")

	button.addEventListener("click", () => {
		editor_state_set(old_state => ({
			notebook: {
				...old_state.notebook,
				process_status: "no_process",
			},
		})).then(() => {
			window.requestAnimationFrame(() => {
				document.querySelector("#process_status a").click()
			})
		})
	})
	button.innerText = "Restart notebook"

	return button
</script>
""";

# ╔═╡ Cell order:
# ╟─af6ae00d-f032-4743-878b-e575466b6e84
# ╠═a70cef7d-2a2f-4155-bdf3-fec9df94c63f
# ╠═87d66a6f-b40d-4d76-afc6-4b4f086e80a4
# ╠═087614cc-a4f6-4867-a6bd-201d1f2a7fc2
# ╠═c7a51cb9-790a-4818-95d7-c32f2252e8b3
# ╠═d87781e8-af89-41d0-a56d-94744406042d
# ╠═4ef3d13c-2f57-4575-bc88-c8092ae6910f
# ╠═4d78d81c-f1a3-4386-adce-2a7fa1e3eaef
# ╟─3a99a9bb-7862-41f8-a535-d3c580da6909
# ╠═b24b7e23-61ea-41fc-a345-286e904c042b
# ╠═1bb47749-edde-4bee-be9f-059a7652b354
# ╠═6aae7cca-01f1-45be-9d4b-a3942bd042f3
# ╠═27c9e166-df00-4a30-a6ba-32e9a071136b
# ╠═773e6c0a-08ea-47b2-8924-42759f944c6b
# ╠═ceaff98b-1ba6-43b3-bd2f-0e7fe585f428
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╟─afe4745f-f9f1-4e23-8735-cbec6fb79c41
# ╟─b8fd36a7-d8d1-45f7-b66e-df9132168bfc

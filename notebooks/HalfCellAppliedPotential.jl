### A Pluto.jl notebook ###
# v0.20.21

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
    using JuliaMPBSolver.ICMPBP: ICMPBP, ICMPBData, AppliedPotentialHalfCell, set_molarity!, calc_cmol, calc_c0mol, calc_χ, get_E, get_φ, get_p, get_c0,
        set_κ!, set_φ!
    using DrWatson, PoissonBoltzmannIPAM2025

end

# ╔═╡ af6ae00d-f032-4743-878b-e575466b6e84
md"""
# Half cell with applied potential
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
    const L = 10nm
    const n = 101
end

# ╔═╡ c7a51cb9-790a-4818-95d7-c32f2252e8b3
Xhalf = (10 .^ range(-3, 1, length = 101)) * nm;

# ╔═╡ 4ef3d13c-2f57-4575-bc88-c8092ae6910f
gridhalf = simplexgrid(Xhalf)

# ╔═╡ 3a99a9bb-7862-41f8-a535-d3c580da6909
md"""
All values are given with respect to SI basic units (m, kg, s, V, A)
"""

# ╔═╡ b24b7e23-61ea-41fc-a345-286e904c042b
datavhalf = ICMPBData(χvar = true)

# ╔═╡ 1bb47749-edde-4bee-be9f-059a7652b354
begin
    datavhalf.Escale = 1
    halfcell = AppliedPotentialHalfCell(gridhalf, datavhalf, dielectric_decrement = false, valuetype = Float64)

    halfcelldd = AppliedPotentialHalfCell(gridhalf, datavhalf, dielectric_decrement = true, valuetype = Float64)
end;


# ╔═╡ dd3c4807-3972-4e9c-a44a-3347b065d01c
p3 = plotcells(halfcell, halfcelldd); p3

# ╔═╡ 935f8897-4d68-447b-8316-1a0fe0285d54
@test isa(p3, Figure)

# ╔═╡ d07ac411-7985-4b5f-a88b-8aa4037b7d65
function makecolors(V)
    h = 0.5 / length(V)
    return colors = [ (i * h, 0, 1 - i * h) for i in 1:length(V) ]
end

# ╔═╡ 2cf54b71-d99b-40bd-b9a6-0a7cf919614b
p1 = let
    molarities = [0.001, 0.01, 0.1, 1]
    φ_max = 0.5
    colors = makecolors(molarities)
    pyplot.close()
    clf()
    fig, ax = pyplot.subplots(1, 1)
    fig.set_size_inches(600 / 100, 300 / 100)
    κ = 10
    ax.set_title("Double layer capacitance, κ=$(κ)\n Dashed: witghout dielectric decrement")
    set_κ!(halfcell, κ)
    set_κ!(halfcelldd, κ)
    for i in 1:length(molarities)
        M = molarities[i]
        set_molarity!(halfcell, M)

        volts, dlcaps = ICMPBP.dlcapsweep(halfcell; φ_max)
        plot(volts, dlcaps / (μF / cm^2), color = colors[i], linestyle = "--")

        set_molarity!(halfcelldd, M)
        volts, dlcaps = ICMPBP.dlcapsweep(halfcelldd; φ_max, damp_initial = 0.1)
        plot(volts, dlcaps / (μF / cm^2), label = "M=$(M)", color = colors[i])
    end
    ax.set_xlabel("U/V")
    ax.set_ylabel(L"C_{dl}/(μF/cm^2)")
    ax.legend()
    ax.grid()
    tight_layout()

    !haskey(ENV, "CI") && savefig(draftresultsdir("halfcell_dlcap_M"), dpi = 600)

    gcf()
end; p1

# ╔═╡ 52d5fb3c-072d-4a6e-a392-7e236b4ec933
@test isa(p1, Figure)

# ╔═╡ cc8e8ba6-4cc6-4e8a-bbd6-07b9f61c2ef5
p2 = let
    kappas = [1, 5, 10, 20]
    φ_max = 0.5
    M = 0.1
    colors = makecolors(kappas)
    pyplot.close()
    clf()
    fig, ax = pyplot.subplots(1, 1)
    fig.set_size_inches(600 / 100, 300 / 100)
    ax.set_title("Double layer capacitance, M=$(M)\n Dashed: without dielectric decrement")
    set_molarity!(halfcell, M)
    set_molarity!(halfcelldd, M)
    for i in 1:length(kappas)
        κ = kappas[i]
        set_κ!(halfcell, κ)
        set_κ!(halfcelldd, κ)

        volts, dlcaps = ICMPBP.dlcapsweep(halfcelldd; φ_max, damp_initial = 0.1)
        plot(volts, dlcaps / (μF / cm^2), label = "κ=$(κ)", color = colors[i])

        volts, dlcaps = ICMPBP.dlcapsweep(halfcell; φ_max)
        plot(volts, dlcaps / (μF / cm^2), color = colors[i], linestyle = "--")

    end
    ax.set_xlabel("U/V")
    ax.set_ylabel(L"C_{dl}/(μF/cm^2)")
    ax.legend()
    ax.grid()
    tight_layout()
    !haskey(ENV, "CI") && savefig(draftresultsdir("halfcell_dlcap_kappa"), dpi = 600)

    gcf()
end; p2

# ╔═╡ 48896662-d87f-4606-9118-6471184b4dc7
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
# ╠═4ef3d13c-2f57-4575-bc88-c8092ae6910f
# ╟─3a99a9bb-7862-41f8-a535-d3c580da6909
# ╠═b24b7e23-61ea-41fc-a345-286e904c042b
# ╠═1bb47749-edde-4bee-be9f-059a7652b354
# ╠═2cf54b71-d99b-40bd-b9a6-0a7cf919614b
# ╠═52d5fb3c-072d-4a6e-a392-7e236b4ec933
# ╠═cc8e8ba6-4cc6-4e8a-bbd6-07b9f61c2ef5
# ╠═48896662-d87f-4606-9118-6471184b4dc7
# ╠═dd3c4807-3972-4e9c-a44a-3347b065d01c
# ╠═935f8897-4d68-447b-8316-1a0fe0285d54
# ╠═d07ac411-7985-4b5f-a88b-8aa4037b7d65
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╟─afe4745f-f9f1-4e23-8735-cbec6fb79c41
# ╟─b8fd36a7-d8d1-45f7-b66e-df9132168bfc

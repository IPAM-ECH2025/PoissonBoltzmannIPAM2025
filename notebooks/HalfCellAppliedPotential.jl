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
    using PythonPlot
	using PythonPlot: pyplot
	using LaTeXStrings
    using Colors
    using JuliaMPBSolver.ICMPBP: ICMPBP, ICMPBData, AppliedPotentialHalfCell, set_molarity!, calc_cmol, calc_c0mol, calc_χ
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
	const L=10nm
	const n=101
end

# ╔═╡ c7a51cb9-790a-4818-95d7-c32f2252e8b3
Xhalf=(10 .^range(-3,1, length=101))*nm;

# ╔═╡ 4ef3d13c-2f57-4575-bc88-c8092ae6910f
gridhalf=simplexgrid(Xhalf)

# ╔═╡ b24b7e23-61ea-41fc-a345-286e904c042b
datavhalf = ICMPBData()

# ╔═╡ 1bb47749-edde-4bee-be9f-059a7652b354
sysvhalf=AppliedPotentialHalfCell(gridhalf, datavhalf);

# ╔═╡ 68380bc6-235b-42d1-bdc8-c843e3f91dca
let
    pyplot.close()
	clf()
	fig, ax = pyplot.subplots(3, 1)
    fig.set_size_inches(600/ 100, 600 / 100)
    ax1 = ax[0]
    ax2 = ax[1]
    ax3 = ax[2]
	κ=20
	M=0.1
	datavhalf.κ=[κ,κ]
	set_molarity!(datavhalf, M)
	fig.suptitle("Solutions for κ=$(κ), M=$(M)\nDashed: with dielectric decrement")

	X=Xhalf/nm
	datavhalf.φ=0.0

	datavhalf.χvar=false
	solc=solve(sysvhalf)
	for φ in [0.1,0.2,0.3, 0.4, 0.5]
		datavhalf.φ=φ
		solc=solve(sysvhalf, inival=solc)
	end
    cc = calc_cmol(solc, sysvhalf.sys)
    c0c = calc_c0mol(solc, sysvhalf.sys)
	φc=solc[datavhalf.iφ,:]
	pc=solc[datavhalf.ip,:]
    Ec=solc[datavhalf.iE, :]
	χc=calc_χ(solc, sysvhalf.sys)

	
	datavhalf.χvar=true
	solv=solve(sysvhalf, inival=solc)
	cv = calc_cmol(solv, sysvhalf.sys)
    c0v = calc_c0mol(solv, sysvhalf.sys)
	φv=solv[datavhalf.iφ,:]
	pv=solv[datavhalf.ip,:]
    Ev=solv[datavhalf.iE, :]
	χv=calc_χ(solv, sysvhalf.sys)
	

	
	
    ax1.grid()
    ax1r = ax1.twinx()
    ax1.semilogx(X, φv, color = "green", linestyle="--")
    ax1r.semilogx(X,pv, color = "red", linestyle="--")
    ax1.semilogx(X, φc, color = "green", label = "ϕ/V")
    ax1r.semilogx(X,pc, color = "red",label = L"p/GPa")
    ax1.set_xlabel("z/nm")
    ax1.legend(loc = (0.2, 0.1))
    ax1.set_ylabel("ϕ/V")
    ax1r.set_ylabel("p/GPa")
    ax1r.legend(loc = (0.8, 0.1))
	ax1r.set_ylim((0,0.2))
    ax2.grid()

    ax2.set_xlabel("z/nm")
    ax2.set_ylabel("c/(mol/L)")


    ax2.semilogx(X, cc[1, :], color = "blue", label = L"c^-")
    ax2.semilogx(X, cc[2, :], color = "red",label = L"c^+")
    ax2.semilogx(X, c0c, color = "green", label = L"c_{solvent}")

	ax2.semilogx(X, cv[1, :], color = "blue", linestyle="--")
    ax2.semilogx(X, cv[2, :], color = "red",linestyle="--")
    ax2.semilogx(X, c0v, color = "green", linestyle="--")
    ax2.legend()


    ax3.grid()
    ax3r = ax3.twinx()
    ax3.semilogx(X, Ec/(V/nm), color = "green", label = "E")
    ax3r.semilogx(X, χc, color = "red", label = "χ")
    ax3.semilogx(X, Ev/(V/nm), color = "green", linestyle="--")
    ax3r.semilogx(X, χv, color = "red", linestyle="--")
    ax3.legend(loc = (0.1, 0.1))
    ax3r.legend(loc = (0.8, 0.5))
	  ax3.set_xlabel("z/nm")
    ax3.set_ylabel("E/(V/nm)")
    ax3r.set_ylabel("χ")
	ax3r.set_ylim((0,100))

    tight_layout()
    gcf()
	
end

# ╔═╡ d07ac411-7985-4b5f-a88b-8aa4037b7d65
function makecolors(V)
	h=0.5/length(V)
	colors=[ (i*h,0,1-i*h) for i in 1:length(V) ]
end

# ╔═╡ 2cf54b71-d99b-40bd-b9a6-0a7cf919614b
let
	molarities=[0.001,0.01, 0.1, 1]
	colors=makecolors(molarities)
	pyplot.close()
	clf()
	fig, ax = pyplot.subplots(1, 1)
    fig.set_size_inches(600 / 100, 300 / 100)
	ax.set_title("Double layer capacitance, κ=10\n Dashed: dielectric decrement")
			datavhalf.κ=[10,10]
	for i in 1:length(molarities)
		M=molarities[i]
		datavhalf.χvar=false
		set_molarity!(datavhalf, M)
		volts, dlcaps=ICMPBP.dlcapsweep(sysvhalf)
		plot(volts,dlcaps/(μF/cm^2), label="M=$(M)", color=colors[i])
		datavhalf.χvar=true
		volts, dlcaps=ICMPBP.dlcapsweep(sysvhalf)
		plot(volts,dlcaps/(μF/cm^2), color=colors[i], linestyle="--")
	end
	ax.set_xlabel("U/V")
	ax.set_ylabel(L"C_{dl}/(μF/cm^2)")
	ax.legend()
	ax.grid()
	gcf()
end

# ╔═╡ cc8e8ba6-4cc6-4e8a-bbd6-07b9f61c2ef5
let
	kappas=[0,5,10,20]
	colors=makecolors(kappas)
	pyplot.close()
	clf()
	fig, ax = pyplot.subplots(1, 1)
    fig.set_size_inches(600 / 100, 300 / 100)
	ax.set_title("Double layer capacitance, M=0.1\n Dashed: dielectric decrement")
		set_molarity!(datavhalf, 0.1)
	for i in 1:length(kappas)
		κ=kappas[i]
		datavhalf.κ=[κ,κ]
		datavhalf.χvar=false
		volts, dlcaps=ICMPBP.dlcapsweep(sysvhalf)
		plot(volts,dlcaps/(μF/cm^2), label="κ=$(κ)", color=colors[i])
		datavhalf.χvar=true
		volts, dlcaps=ICMPBP.dlcapsweep(sysvhalf)
		plot(volts,dlcaps/(μF/cm^2), color=colors[i], linestyle="--")
	end
	ax.set_xlabel("U/V")
	ax.set_ylabel(L"C_{dl}/(μF/cm^2)")
	ax.legend()
	ax.grid()
	gcf()
end

# ╔═╡ 8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
html"""<hr>"""

# ╔═╡ 784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
html"""<style>.dont-panic{ display: none }</style>"""

# ╔═╡ afe4745f-f9f1-4e23-8735-cbec6fb79c41
begin
    function floataside(text::Markdown.MD; top = 1)
        uuid = uuid1()
        @htl("""
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

       		""")
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
# ╠═b24b7e23-61ea-41fc-a345-286e904c042b
# ╟─1bb47749-edde-4bee-be9f-059a7652b354
# ╟─2cf54b71-d99b-40bd-b9a6-0a7cf919614b
# ╟─cc8e8ba6-4cc6-4e8a-bbd6-07b9f61c2ef5
# ╟─68380bc6-235b-42d1-bdc8-c843e3f91dca
# ╟─d07ac411-7985-4b5f-a88b-8aa4037b7d65
# ╟─8af12f1c-d35b-4cc9-8185-1bb5adbb69e8
# ╟─784b4c3e-bb2a-4940-a83a-ed5e5898dfd4
# ╟─afe4745f-f9f1-4e23-8735-cbec6fb79c41
# ╟─b8fd36a7-d8d1-45f7-b66e-df9132168bfc

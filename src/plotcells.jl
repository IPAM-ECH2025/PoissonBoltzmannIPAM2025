surfcharge(n) = n * ph"e" / ufac"nm^2"

function makecolors(V)
    h = 0.5 / length(V)
    return colors = [ (i * h, 0, 1 - i * h) for i in 1:length(V) ]
end

function plotcells(
        cell::AbstractSymmetricCell, celldd::AbstractSymmetricCell;
        figname = "symmcell",
        κ = 10, M = 1, q = 1.0,
        title = latexstring(""" \\text{Solutions for } \\kappa=$(κ), M_{avg}=$(M)mol/L, q=$(q)e/nm^2""") *
            "\n" *
            latexstring("""\\text{Dashed: without dielectric decrement} """)

    )
    pyplot.close()
    clf()
    fig, ax = pyplot.subplots(3, 1)
    fig.set_size_inches(600 / 100, 600 / 100)
    ax1 = ax[0]
    ax2 = ax[1]
    ax3 = ax[2]
    set_κ!(cell, κ)
    set_κ!(celldd, κ)
    set_molarity!(cell, M)
    set_molarity!(celldd, M)
    fig.suptitle(title)
    grid = cell.sys.grid
    X = grid[Coordinates][1, :] / ufac"nm"
    solc = unknowns(cell)
    for _q in range(0, q, length = 5)
        set_q!(cell, surfcharge(_q))
        solc = solve(cell, inival = solc)
    end

    cc = calc_cmol(solc, cell)
    c0c = calc_c0mol(solc, cell)
    φc = get_φ(solc, cell)
    pc = get_p(solc, cell) / ufac"GPa"
    Ec = get_E(solc, cell) .|> abs
    χc = calc_χ(solc, cell)


    solv = unknowns(celldd)
    for _q in range(0, q, length = 10)
        set_q!(celldd, surfcharge(_q))
        solv = solve(celldd, inival = solv, damp_initial = 0.1)
    end
    cv = calc_cmol(solv, celldd)
    c0v = calc_c0mol(solv, celldd)
    φv = get_φ(solv, celldd)
    pv = get_p(solv, celldd) / ufac"GPa"
    Ev = get_E(solv, celldd) .|> abs
    χv = calc_χ(solv, celldd)

    ax1r = ax1.twinx()
    ax1.plot(X, φc, color = "olive", linestyle = "--")
    ax1r.plot(X, pc, color = "darkorange", linestyle = "--")
    ax1.plot(X, φv, color = "olive", label = "ϕ")
    ax1r.plot(X, pv, color = "darkorange", label = "p")
    ax1.set_xlabel("z/nm")
    ax1.legend(loc = (0.2, 0.1))
    ax1.set_ylabel("ϕ/V")
    φ_max = ceil(maximum(φv) * 10) / 10
    ax1r.set_ylabel("p/GPa")
    ax1r.legend(loc = (0.7, 0.1))
    ax1.set_ylim((-φ_max, φ_max))
    ax1.set_yticks(range(-φ_max, φ_max, length = 5))

    pmax = ceil(20 * maximum(pc)) / 10
    ax1r.set_ylim((0, pmax))
    ax1r.set_yticks(range(0, pmax, length = 5))
    ax1.grid()


    ax2.grid()

    ax2.set_xlabel("z/nm")
    ax2.set_ylabel("c/(mol/L)")
    ax2r = ax2.twinx()
    ax2r.set_ylabel(L"c_{solvent}/(mol/L)")


    ax2.plot(X, cv[1, :], color = "blue", label = L"c^-")
    ax2.plot(X, cv[2, :], color = "red", label = L"c^+")
    ax2r.plot(X, c0v, color = "green", label = L"c_{solvent}")

    ax2.plot(X, cc[1, :], color = "blue", linestyle = "--")
    ax2.plot(X, cc[2, :], color = "red", linestyle = "--")
    ax2r.plot(X, c0c, color = "green", linestyle = "--")
    ax2.legend(loc = (0.15, 0.6))
    ax2r.legend(loc = (0.7, 0.7))
    ax2.set_ylim((0, 5))
    ax2.set_yticks(range(0, 5, length = 5))
    ax2r.set_ylim((0, 60))
    ax2r.set_yticks(range(0, 60, length = 5))

    ax3.grid()
    ax3r = ax3.twinx()
    ax3.plot(X, Ev / ufac"V / nm", color = "darkviolet", label = "|E|")
    ax3r.plot(X, χv, color = "steelblue", label = "χ")
    ax3.plot(X, Ec / ufac"V / nm", color = "darkviolet", linestyle = "--")
    ax3.legend(loc = (0.2, 0.2))
    ax3r.plot(X, χc, color = "steelblue", linestyle = "--")
    ax3r.legend(loc = (0.65, 0.2))
    ax3.set_xlabel("z/nm")
    Emax = ceil(maximum(Ev / ufac"V / nm"))
    ax3.set_ylabel("|E|/(V/nm)")
    ax3.set_ylim((0, Emax))
    ax3.set_yticks(range(0, Emax, length = 5))
    ax3r.set_ylabel("χ")
    ax3r.set_ylim((0, 100))
    ax3r.set_yticks(range(0, 100, length = 5))

    tight_layout()
    !haskey(ENV, "CI") && savefig(draftresultsdir(figname), dpi = 600)
    return gcf()

end


function plotcells(
        cell::AbstractHalfCell, celldd::AbstractHalfCell;
        figname = "halfcell",
        κ = 10, M = 1, Δφ = 0.5,
        title = latexstring(""" \\text{Solutions for } κ=$(κ), M=$(M)mol/L, \\Delta \\varphi=$(Δφ)V""") *
            "\n" *
            latexstring("""\\text{Dashed: without dielectric decrement} """)
    )
    pyplot.close()
    clf()
    fig, ax = pyplot.subplots(3, 1)
    fig.set_size_inches(600 / 100, 600 / 100)
    ax1 = ax[0]
    ax2 = ax[1]
    ax3 = ax[2]
    set_κ!(cell, κ)
    set_κ!(celldd, κ)
    set_molarity!(cell, M)
    set_molarity!(celldd, M)
    fig.suptitle(title)
    grid = cell.sys.grid
    X = grid[Coordinates][1, :] / ufac"nm"
    solc = unknowns(cell)
    for _Δφ in range(0, Δφ, length = 5)
        set_φ!(cell, _Δφ)
        solc = solve(cell, inival = solc)
    end

    cc = calc_cmol(solc, cell)
    c0c = calc_c0mol(solc, cell)
    φc = get_φ(solc, cell)
    pc = get_p(solc, cell) / ufac"GPa"
    Ec = get_E(solc, cell) .|> abs
    χc = calc_χ(solc, cell)


    solv = unknowns(celldd)
    set_φ!(celldd, 0)
    solv = solve(celldd, inival = solv, damp_initial = 0.5)
    h = 0.001
    v = 0.0
    while v < Δφ
        try
            vv = min(v + h, Δφ)
            set_φ!(celldd, vv)
            solv = solve(celldd, inival = solv, damp_initial = 0.5)
            v = vv
            h = h * 1.2
        catch e
            h = h / 2
            if h < 1.0e-5
                rethrow(e)
            end
        end
    end
    #set_φ!(celldd, Δφ)
    #solv = solve(celldd, inival = solc, damp_initial = 0.1)

    cv = calc_cmol(solv, celldd)
    c0v = calc_c0mol(solv, celldd)
    φv = get_φ(solv, celldd)
    pv = get_p(solv, celldd) / ufac"GPa"
    Ev = get_E(solv, celldd) .|> abs
    χv = calc_χ(solv, celldd)

    ax1r = ax1.twinx()
    ax1.plot(X, φc, color = "olive", linestyle = "--")
    ax1r.plot(X, pc, color = "darkorange", linestyle = "--")
    ax1.plot(X, φv, color = "olive", label = "ϕ")
    ax1r.plot(X, pv, color = "darkorange", label = L"p")
    ax1.set_xlabel("z/nm")
    ax1.legend(loc = (0.2, 0.7))
    ax1.set_ylabel("ϕ/V")
    φ_max = ceil(maximum(φv) * 10) / 10
    ax1r.set_ylabel("p/GPa")
    ax1r.legend(loc = (0.7, 0.7))
    ax1.set_ylim((0, φ_max))
    ax1.set_yticks(range(0, φ_max, length = 5))

    pmax = ceil(20 * maximum(pc)) / 10
    ax1r.set_ylim((0, pmax))
    ax1r.set_yticks(range(0, pmax, length = 5))
    ax1.grid()


    ax2.grid()

    ax2.set_xlabel("z/nm")
    ax2.set_ylabel("c/(mol/L)")
    ax2r = ax2.twinx()
    ax2r.set_ylabel(L"c_{solvent}/(mol/L)")


    ax2.plot(X, cv[1, :], color = "blue", label = L"c^-")
    ax2.plot(X, cv[2, :], color = "red", label = L"c^+")
    ax2r.plot(X, c0v, color = "green", label = L"c_{solvent}")

    ax2.plot(X, cc[1, :], color = "blue", linestyle = "--")
    ax2.plot(X, cc[2, :], color = "red", linestyle = "--")
    ax2r.plot(X, c0c, color = "green", linestyle = "--")
    ax2.legend(loc = (0.2, 0.6))
    ax2r.legend()
    ax2.set_ylim((0, 5))
    ax2.set_yticks(range(0, 5, length = 5))
    ax2r.set_ylim((0, 60))
    ax2r.set_yticks(range(0, 60, length = 5))

    ax3.grid()
    ax3r = ax3.twinx()
    ax3.plot(X, Ev / ufac"V / nm", color = "darkviolet", label = "E")
    ax3r.plot(X, χv, color = "steelblue", label = "χ")
    ax3.plot(X, Ec / ufac"V / nm", color = "darkviolet", linestyle = "--")
    ax3.legend(loc = (0.2, 0.1))
    ax3r.plot(X, χc, color = "steelblue", linestyle = "--")
    ax3r.legend(loc = (0.8, 0.1))
    ax3.set_xlabel("z/nm")
    Emax = ceil(maximum(Ev / ufac"V / nm"))
    ax3.set_ylabel("|E|/(V/nm)")
    ax3.set_ylim((0, Emax))
    ax3.set_yticks(range(0, Emax, length = 5))
    ax3r.set_ylabel("χ")
    ax3r.set_ylim((0, 100))
    ax3r.set_yticks(range(0, 100, length = 5))

    tight_layout()
    !haskey(ENV, "CI") && savefig(draftresultsdir(figname), dpi = 600)
    return gcf()

end

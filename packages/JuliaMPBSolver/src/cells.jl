abstract type AbstractMPBCell end

abstract type AbstractHalfCell <: AbstractMPBCell end
abstract type AbstractSymmetricCell <: AbstractMPBCell end


function VoronoiFVM.unknowns(cell::AbstractMPBCell)
    sys = cell.sys
    data = sys.physics.data
    (; i0, iφ, ip, iE, coffset, N) = data
    u = unknowns(sys; inival = 0)
    i3 = 0
    if data.conserveions
        i3 = sys.grid[BFaceNodes][3][1]
    end
    for α in 1:N
        u[α, :] .= 0.1
        if data.conserveions
            u[coffset + α, i3] = data.n_E[α] / data.cscale
        end
    end

    u[i0, :] .= 1 - N * 0.1
    return u
end

mpbdata(cell) = cell.sys.physics.data

calc_cmol(sol, cell::AbstractMPBCell) = calc_cmol(sol, cell.sys)
calc_c0mol(sol, cell::AbstractMPBCell) = calc_c0mol(sol, cell.sys)
calc_χ(sol, cell::AbstractMPBCell) = calc_χ(sol, cell.sys)
get_E(sol, cell::AbstractMPBCell) = sol[mpbdata(cell).iE, :] * mpbdata(cell).Escale
get_φ(sol, cell::AbstractMPBCell) = sol[mpbdata(cell).iφ, :]
get_p(sol, cell::AbstractMPBCell) = sol[mpbdata(cell).ip, :] * mpbdata(cell).pscale
get_c0(sol, cell::AbstractMPBCell) = sol[mpbdata(cell).i0, :]
set_κ!(cell::AbstractMPBCell, κ::Number) = mpbdata(cell).κ = [κ, κ]
set_molarity!(cell::AbstractMPBCell, M) = set_molarity!(mpbdata(cell), M)
set_φ!(cell::AbstractMPBCell, φ::Number) = mpbdata(cell).φ = φ

function SciMLBase.solve(cell::AbstractMPBCell; inival = unknowns(cell), verbose = "", damp_initial = 0.1, kwargs...)
    sys = cell.sys
    return solve(sys; inival, damp_initial, verbose, kwargs...)
end

struct AppliedPotentialHalfCell <: AbstractHalfCell
    sys::VoronoiFVM.System
end


function halfcell_applied_potential_bcondition!(y, u, bnode, data)
    (; iφ, ip) = data
    boundary_dirichlet!(y, u, bnode, species = iφ, region = 1, value = data.φ)
    boundary_dirichlet!(y, u, bnode, species = iφ, region = 2, value = 0)
    boundary_dirichlet!(y, u, bnode, species = ip, region = 2, value = 0)
    return nothing
end

function AppliedPotentialHalfCell(grid, data; dielectric_decrement = false, valuetype = Float64)
    data = deepcopy(data)
    data.nv = ones(num_nodes(grid)) # help to satisfy sparsity detector
    data.conserveions = false
    data.χvar = dielectric_decrement
    data.χvarext = dielectric_decrement

    sys = VoronoiFVM.System(
        grid;
        data = data,
        flux = poisson_and_p_flux!,
        reaction = reaction!,
        bcondition = halfcell_applied_potential_bcondition!,
        generic = ionconservation!,
        unknown_storage = :dense,
        valuetype
    )

    # Enable species for all fields
    for i in 1:data.N
        enable_species!(sys, i, [1])
    end

    enable_species!(sys, data.i0, [1])
    enable_species!(sys, data.iφ, [1])
    enable_species!(sys, data.ip, [1])
    enable_species!(sys, data.iE, [1])
    data.nv = nodevolumes(sys)
    return AppliedPotentialHalfCell(sys)
end

function dlcapsweep(cell::AppliedPotentialHalfCell; φ_max = 1.0, δφ = 1.0e-3, steps = 101, damp_initial = 1, kwargs...)
    sys = cell.sys
    data = sys.physics.data
    apply_voltage!(data, 0)
    sol0 = solve(cell; damp_initial)
    volts = zeros(0)
    dlcaps = zeros(0)
    for dir in [-1, 1]
        sol = sol0
        for φ in range(0, φ_max, length = steps)

            apply_voltage!(data, dir * φ)
            sol = solve(cell; inival = sol, damp_initial, kwargs...)

            apply_voltage!(data, dir * (φ + δφ))
            solδ = solve(cell; inival = sol, damp_initial, kwargs...)

            Q = calc_spacecharge(sys, sol)
            Qδ = calc_spacecharge(sys, solδ)
            cdl = (Qδ - Q) / (dir * δφ)
            push!(volts, dir * φ)
            push!(dlcaps, cdl)
        end
        if dir == -1
            volts = reverse(volts)[1:(end - 1)]
            dlcaps = reverse(dlcaps)[1:(end - 1)]
        end
    end
    return volts, dlcaps
end


struct AppliedPotentialSymmetricCell <: AbstractMPBCell
    sys::VoronoiFVM.System
end

struct SurfaceChargedHalfCell <: AbstractMPBCell
    sys::VoronoiFVM.System
end
struct SurfaceChargedSymmetricCell <: AbstractMPBCell
    sys::VoronoiFVM.System
end

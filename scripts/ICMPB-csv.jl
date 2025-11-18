using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))
using Revise
using PlutoUI, HypertextLiteral, UUIDs
using LinearAlgebra
using Interpolations
using VoronoiFVM, ExtendableGrids
using LaTeXStrings
using LessUnitful, Unitful
using PreallocationTools
using LaTeXStrings
using DoubleFloats
using ForwardDiff
using PythonPlot
using NLsolve
using DelimitedFiles

using DrWatson

using PoissonBoltzmannIPAM2025
using JuliaMPBSolver

const χ_S = 78.49 - 1
const F = ph"N_A" * ph"e"
const K = ufac"K"
const nm = ufac"nm"
const m = ufac"m"
const dm = ufac"dm"
const V = ufac"V"
const mol = ufac"mol"
const T = (273.15 + 25) * ufac"K"
const RT = ph"R" * T
const c̄ = 55.508mol / dm^3 # summary molar concentration
const ε_0 = ph"ε_0"

const f_mod = true # model choice: 0: Boltzmann, 1: Bikerman
const nref = 4 # grid refinement level
const L = 2.278592867nm # computational domain size
const n_e = 1 # number of electrons/nm^2 at interfaces, defines q
const M_avg = 2.0 # average molarity
const E_0 = 10V / nm # decrement parameter
const a = 5.0 / E_0^2 # decrement parameter in χ(E)
const z = [-1, 1] # species charge numbers

const N = length(z)

const c_avg = fill(M_avg * mol / dm^3, N)

surfcharge(n) = n * ph"e" / ufac"nm^2"

const q = surfcharge(n_e) # surface charge

const c0_avg = c̄ - sum(c_avg) # solvent bulk molar concentration

const l_debye = sqrt((1 + χ_S) * ε_0 * RT / (F^2 * c_avg[1])) |> u"nm"

begin
  X = range(0, L, length = 20 * 2^nref + 1)
  grid = simplexgrid(X)
  bfacemask!(grid, [L / 2], [L / 2], 3, tol = 1.0e-10 * nm)
end

const i3 = grid[BFaceNodes][1, 3] # Index  of grid midpoint

Base.@kwdef mutable struct PBData
  c_ref::Vector{Float64} = c_avg
  c_avg::Vector{Float64} = c_avg
  c̄::Float64 = c̄
  c0_ref::Float64 = c̄ - sum(c_ref)
  q::Float64 = q
  z::Vector{Float64} = z
  N::Int = length(z)
  a::Float64 = a
  E_0::Float64 = E_0
  f_mod::Bool = f_mod
  F::Float64 = F
  RT::Float64 = RT
  ε_0::Float64 = ε_0
  χ_S::Float64 = χ_S
  iϕ::Int = 1
  cache = DiffCache(ones(Float64, length(z)))
end

M3_avg = 2.0
n3_e = 1

function molfractions!(y, ϕ, c_ref, c0_ref, data)
  (; z, N, f_mod, F, RT) = data
  for ic in 1:N
    y[ic] = exp(-z[ic] * ϕ * F / RT) * c_ref[ic] / c0_ref
  end
  denom = 1.0 / (one(ϕ) + f_mod * sum(y))
  for ic in 1:N
    y[ic] = y[ic] * denom
  end
  return nothing
end

function concentrations(sol, data; c_ref = data.c_ref)
  (; c̄, iϕ, N) = data
  n = size(sol, 2)
  c = zeros(N, n)
  y = zeros(N)
  c0_ref = data.c̄ - sum(c_ref)
  for iz in 1:n
    molfractions!(y, sol[data.iϕ, iz], c_ref, c0_ref, data)
    for ic in 1:N
      c[ic, iz] = c̄ * y[ic]
    end
  end
  return c
end

function avgconcentrations(sol, sys; data = data(sys), c_ref = data.c_ref)
  c = concentrations(sol, data; c_ref)
  (; N = data)
  cavg = zeros(N)
  nv = nodevolumes(sys)
  L = sum(nv)
  for ic in 1:N
    cavg[ic] = nv ⋅ c[ic, :] / L
  end
  return cavg
end

function flux!(y, u, edge, data)
  (; χ_S, a, ε_0, iϕ) = data
  eins = one(eltype(u))
  h = edgelength(edge)
  E = (u[iϕ, 1] - u[iϕ, 2]) / h
  χ = χ_S / sqrt(eins + a * E^2)
  ε = (eins + χ) * ε_0
  y[iϕ] = ε * ((u[iϕ, 1] - u[iϕ, 2]))
  return nothing
end

function relpermittivity(sol, data; grid)
  (; χ_S, a, iϕ) = data
  X = grid[Coordinates][1, :]
  return χ_S ./ (
    a *
    ((sol[iϕ, 2:end] - sol[iϕ, 1:(end-1)]) ./ (X[2:end] - X[1:(end-1)])) .^ 2 .+
    1
  ) .+ 1
end

function spacecharge!(y, ϕ, c_ref, c0_ref, data)
  (; N, F, c_ref, c0_ref, c̄, z) = data
  molfractions!(y, ϕ, c_ref, c0_ref, data)
  sumyz = zero(ϕ)
  for i in 1:N
    sumyz += z[i] * y[i]
  end
  return F * c̄ * sumyz
end

function bcondition!(y, u, bnode, data)
  boundary_neumann!(y, u, bnode, species = data.iϕ, region = 2, value = -data.q)
  boundary_neumann!(y, u, bnode, species = data.iϕ, region = 1, value = data.q)
  return nothing
end

function extcref(cref0, data)
  (; z, N) = data
  return push!(copy(cref0), -z[1:(N-1)] ⋅ cref0 / z[N])
end

function clonedata(data0, c_ref)
  data = deepcopy(data0)
  data.c_ref .= c_ref
  data.c0_ref = c̄ - sum(data.c_ref)
  return data
end

function ICMPBSystem(; data = PBData(), generic = VoronoiFVM.nofunc)
  sys = VoronoiFVM.System(
    grid;
    data,
    generic,
    flux = flux!,
    bcondition = bcondition!,
    unknown_storage = :sparse,
  )
  enable_species!(sys, data.iϕ, [1])
  for ic in 1:(data.N-1)
    enable_boundary_species!(sys, data.iϕ + ic, [3])
  end
  return sys
end

sys3_0 = ICMPBSystem();

const nv = nodevolumes(sys3_0);

const idx = unknown_indices(unknowns(sys3_0));

function xreaction!(f, u, sys, data)
  (; cache, iϕ, N, z, F, c_avg, c̄) = data
  y = get_tmp(cache, u)
  c_ref = [u[idx[ic+iϕ, i3]] for ic in 1:(N-1)]
  push!(c_ref, -z[1:(N-1)] ⋅ c_ref / z[N])
  c0_ref = c̄ - sum(c_ref)
  L = sum(nv)

  for ic in 1:(N-1)
    f[idx[ic+iϕ, i3]] = 0
  end

  for iv in 1:length(nv)
    molfractions!(y, u[idx[iϕ, iv]], c_ref, c0_ref, data)
    f[idx[iϕ, iv]] = 0
    for ic in 1:N
      f[idx[iϕ, iv]] -= y[ic] * c̄ * nv[iv] * z[ic] * F
    end
    for ic in 1:(N-1)
      f[idx[ic+iϕ, i3]] += y[ic] * c̄ * nv[iv]
    end
  end
  for ic in 1:(N-1)
    f[idx[ic+iϕ, i3]] = f[idx[ic+iϕ, i3]] - c_avg[ic] * L
  end
  return nothing
end

Z = grid[Coordinates] / nm

indata = (M = M3_avg, q = n3_e, L = L / nm, n = length(Z))

savename(indata)

c3_avg = fill(M3_avg * mol / dm^3, 2)

data3 = PBData(c_avg = c3_avg, c_ref = 0.5 * c3_avg, q = surfcharge(n3_e))

sys3 = ICMPBSystem(data = data3, generic = xreaction!)

state3 = VoronoiFVM.SystemState(sys3; data = data3)

begin
  inival3 = unknowns(sys3, inival = 0.0)
  inival3[2:N, i3] .= c3_avg[1:(N-1)] / 2
end;

sol3 = solve!(state3; inival = inival3, verbose = "n", damp_initial = 0.5)

conc =
  concentrations(sol3, data3, c_ref = extcref(sol3[2:N, i3], data3)) /
  (mol / dm^3)

c0 = c̄ / (mol / dm^3) .- sum(conc, dims = 1)

result = vcat(Z, c0, conc)'

reference_data_filename = resultsdir(savename("icmpb", indata, "hdf5"))
print(reference_data_filename)

if !isfile(reference_data_filename)
  throw(ErrorException("Reference data not found"))
end

(x_ref, c_solvent_ref) =
  JuliaMPBSolver.DataOut.read_hdf5_data(reference_data_filename, "c_solvent")
(x_ref, c_anion_ref) =
  JuliaMPBSolver.DataOut.read_hdf5_data(reference_data_filename, "c_anion")
(x_ref, c_cation_ref) =
  JuliaMPBSolver.DataOut.read_hdf5_data(reference_data_filename, "c_cation")

same_x = x_ref == JuliaMPBSolver.Grid.get_coordinates(grid)
if !same_x
  throw(ErrorException("Grid coordinates mismatch"))
end
same_solvent_concentration = c_solvent_ref == c0
if !same_solvent_concentration
  throw(ErrorException("Solvent concentration mismatch"))
end
same_anion_concentration = c_anion_ref == conc[1, :]
if !same_anion_concentration
  throw(ErrorException("Anion concentration mismatch"))
end
same_cation_concentration = c_cation_ref == conc[2, :]
if !same_cation_concentration
  throw(ErrorException("Cation concentration mismatch"))
end
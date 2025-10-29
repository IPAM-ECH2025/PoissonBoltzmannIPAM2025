module Grid

using ExtendableGrids
using GridVisualize
using CairoMakie

struct GeometricGrid
  domain_size::Float32
  refinement::Int32
  hmin::Float32
  hmax::Float32
  use_offset::Bool
end

GeometricGrid() = GeometricGrid(10.0f0, 1, 0.1f0, 1.0f0, false)

function create_grid(grid::GeometricGrid)
  # Compute a minimum and maximum cell size from the provided values and
  # refinement. This is a little weird to me because the GeometricGrid
  # hmin and hmax should be obeyed here (logically). However, then
  # the refinement is of little use. 
  # TODO: IDK Fix this logic
  local_hmin = grid.hmin / 2.0^grid.refinement
  local_hmax = grid.hmax / 2.0^grid.refinement

  # Create a little offset 
  offset = 0.0f0
  if grid.use_offset
    offset = 1.0e-3 * grid.domain_size
  end

  # Create to geometric spacings and glue them together. This way
  # the coarsest cells are in the middle
  x_left = geomspace(offset, grid.domain_size / 2.0, local_hmin, local_hmax)
  x_right =
    geomspace(grid.domain_size / 2.0, grid.domain_size, local_hmax, local_hmin)
  x = glue(x_left, x_right)

  return simplexgrid(x)
end

function plot_grid(directory, name = "grid.png"; plotter = CairoMakie)
  # TODO: It seems to plot a window somewhere in here. Need to disable that.
  # TODO: Accept a grid for input so this isn't just a weird debugging function.
  if !isnothing(plotter)
    grid = create_grid(GeometricGrid())
    fig = gridplot(grid; Plotter = plotter, resolution = (500, 500))
    plotter.save(joinpath(directory, name), fig)
  end
  return nothing
end

end
module DataOut

using HDF5
using ExtendableGrids

using ..Grid

function write_data!(filename::String, grid::ExtendableGrid, solution)
  @assert dim_grid(grid) == 1

  file_stream = h5open(filename, "cw")

  grid_coordinates = Grid.get_coordinates(grid)

  @assert length(grid_coordinates) == length(solution)

  file_stream["x"] = grid_coordinates
  file_stream["solution"] = solution

  close(file_stream)
  return nothing
end

function read_data(filename::String)
  file_stream = h5open(filename, "r")

  grid_coordinates = file_stream["x"][:]
  solution = file_stream["solution"][:]

  close(file_stream)

  return grid_coordinates, solution
end

export write_data!

end

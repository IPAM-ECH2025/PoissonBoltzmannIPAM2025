module JuliaMPBSolver

include("units.jl")
using .Units

include("parameters.jl")
using .Parameters

include("grid.jl")
using .Grid

include("postprocess.jl")
using .Postprocess

include("equations.jl")
using .Equations


module ICMPBP
    using LessUnitful
    using ExtendableGrids
    using VoronoiFVM
    using LinearAlgebra
    using SciMLBase
    include("icmbp-p.jl")
    include("cells.jl")
end

include("api.jl")

export mpbpsolve
export icmpbpsolve
end

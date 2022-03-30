module MultilayerCanopy

import FastGaussQuadrature: gausslegendre
using StaticArrays
import LinearAlgebra: UpperTriangular, LowerTriangular
import SpecialFunctions: beta
import QuadGK: quadgk
import HCubature: hcubature

# General
include("Parameters.jl")
include("Integrators.jl")
# Photosynthesis
include("LeafModel.jl")
include("Micrometeo.jl")
include("LeafAngles.jl")
include("CanopyModel.jl")

# Write your package code here.

end

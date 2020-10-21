"""A package for graph computations related to graph Laplacians

Graphs are represented by sparse adjacency matrices, etc.
"""
module Laplacians


  using LinearAlgebra, SparseArrays, Statistics

  include("approxCholTypes.jl")
  include("approxChol.jl")
  include("graphOps.jl")
  include("graphUtils.jl")

  export adj, lap, extendMatrix, approxchol_lap_pc, LDLsolver!, LDLinv

end # module Laplacians.jl

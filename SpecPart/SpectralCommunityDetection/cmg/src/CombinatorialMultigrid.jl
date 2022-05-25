module CombinatorialMultigrid
    using SparseArrays
    using LinearAlgebra
    using LDLFactorizations
    using Laplacians

    include("cmgAlg.jl")
    export cmg_preconditioner_adj, cmg_preconditioner_lap, lPreconditioner
end

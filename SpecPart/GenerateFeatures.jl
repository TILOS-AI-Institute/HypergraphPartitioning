module HypergraphPartitionFeatures

using Shuffle
using JuMP
using Cbc
using LinearAlgebra
using DataStructures
using SparseArrays
using Random
using Statistics
using LightGraphs
using SimpleWeightedGraphs
using DataFrames
using CSV

include("PartitionStructures.jl")
include("SpectralCommunityDetection/SpectralCommunities.jl")
include("GenerateHyperedgesHash.jl")
include("HypergraphToPairList.jl")
include("StructureDetection.jl")
include("HypergraphToIncidence.jl")
include("ReadFiles.jl")
include("PartitionILP.jl")
include("Coarsen.jl")
include("Uncoarsen.jl")
include("GenerateCoarseHypergraph.jl")
include("ExportHypergraph.jl")
include("EstimateClusteringMistakes.jl")
include("InitialPartitions.jl")
include("FM.jl")

function GenerateFeatures(;hg::String = "", nparts::Int = 2, nruns::Int = 10, ub::Int = 5, seed::Int = 0)
    hypergraph_file = hg
    num_parts = nparts
    num_runs = nruns
    ub_factor = ub
    Random.seed!(seed)

    line_log = repeat("=", 80)
    @info "$line_log"
    @info "STARTING FEATURE GENERATING ENGINE"
    @info "$line_log"

    t_elapsed_io = @elapsed begin
        (hedges, eptr, vwts, n, e) = ReadHmetisBenchmark(hypergraph_file)
        hwts = ones(Int, e)
        hypergraph = Hypergraph(n, e, hedges, eptr, vwts, hwts)
        hsizes = hypergraph.eptr[2:end] - hypergraph.eptr[1:end-1]
        hsize_max = maximum(hsizes)
        incidence_struct = HypergraphToIncidence(hypergraph)
        adj_mat = HypergraphToGraph(hypergraph, 1)
    end

    max_capacity = Int(ceil(sum(vwts) * (50+ub_factor)/100))
    min_capacity = sum(vwts) - max_capacity
    capacities = [max_capacity, max_capacity]
    fixed_part = -ones(Int, n)

    @info "TOTAL VERTICES: $n"
    @info "TOTAL HYPEREDGES: $e"
    @info "SIZE OF LARGEST HYPEREDGE: $hsize_max"
    @info "MAX CAPACITY CONSTRAINT: $max_capacity"
    @info "MIN CAPACITY CONSTRAINT: $min_capacity"
    @info "$line_log"

    degree_vector = incidence_struct.eptr[2:end] - incidence_struct.eptr[1:end-1]
    p_0 = findall(x-> x == 0, fixed_part)
    p_1 = findall(x-> x == 1, fixed_part)
    fixed_vertices = Pindex(p_0, p_1)
    eigenvectors = GenEigenVecs(hypergraph, adj_mat[1], fixed_vertices, false, 2, 50)
    hyperstretch_hyperedges = GenerateHyperstretch(hypergraph, eigenvectors)
    incident_hyperstretch_vector = zeros(hypergraph.n, 2)

    for i in 1:hypergraph.n
        start_idx = incidence_struct.eptr[i]
        end_idx = incidence_struct.eptr[i+1]-1
        hdges = incidence_struct.hedges[start_idx:end_idx]
        hdges_hyperstretch = hyperstretch_hyperedges[hdges,:]
        total_incident_hyperstretch = sum(hdges_hyperstretch, dims=1)
        incident_hyperstretch_vector[i, :] = total_incident_hyperstretch
    end

    pvec = zeros(Int, hypergraph.n)

    f = open(hypergraph_file*".part.2")
    i = 0

    for ln in eachline(f)
        i += 1
        pvec[i] = parse(Int, ln)
    end

    close(f)

    f = open(hypergraph_file*".features.dat", "w")

    df = DataFrame([:Eig1 => eigenvectors[:, 1], :Eig2 => eigenvectors[:, 2], :Stretch1 => incident_hyperstretch_vector[:, 1], :Stretch2 => incident_hyperstretch_vector[:, 2], :Degree => degree_vector, :Partition => pvec])
    CSV.write(hypergraph_file*".features.dat", df)

    close(f)

    return (hyperstretch_hyperedges, degree_vector, n)
end

end #module ends
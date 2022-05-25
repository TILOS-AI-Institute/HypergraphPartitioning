module ExactRefiner

using SparseArrays, JuMP, Cbc, GLPK, Statistics

struct Hypergraph
	n::Int
	e::Int
	hedges::Vector{Int}
	eptr::Vector{Int}
	w_::Vector{Int}
	weighted::Bool
end

struct Incidence
	n::Int64
	e::Int64
	hedges::Vector{Int64}
	eptr::Vector{Int64}
	d::Vector{Int64}
end

include("Hypergraph2Incidence.jl")
include("ReadFiles.jl")
include("CutProfile.jl")
include("PartitionILP.jl")
include("Clustering.jl")
include("CombinePartitions.jl")

function xRefine(fname::String, pname::String, nsolns::Int, UBFactor::Int)
    hE = 0
    hedges_picked = Int[]
    bins = Int[]
    hedges = Int[]
    eptr = Int[]
    vwts = Int[]

    (bins, hedges, eptr, vwts, n, hE, ~) = readPartitions(fname, pname)
    vwts = ones(Int, n)
    H = Hypergraph(n, hE, hedges, eptr, ones(Int, hE), false)
    B = hypergraph2incidence(H)
    (cutlist, pids) = prunePartitions(nsolns, bins, H, B)
    (bins_best, fixed) = processBestPartitions(bins, pids)
    #fixed = -ones(Int, n)
    (H_c, B_c, fixed_c, vwts_c, cc) = clusterInstance(bins_best, H, B, vwts, fixed)
    max_capacity = Int(ceil(sum(vwts) * (50+UBFactor)/100))
    capacities = [max_capacity, max_capacity]
    bins_ilp = partitionILP(H_c, vwts_c, capacities, fixed_c)
    hmetis_cut =  cutlist[1][1]
    xRefine_cut = findCutsize(bins_ilp, H_c, B_c)

    return (hmetis_cut, xRefine_cut, H_c.n, H_c.e)
end

end
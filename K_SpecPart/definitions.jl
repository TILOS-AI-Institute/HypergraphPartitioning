struct __hypergraph__
    num_vertices::Int
    num_hyperedges::Int
    eptr::Vector{Int}
    eind::Vector{Int}
    vptr::Vector{Int}
    vind::Vector{Int}
    fixed::Vector{Int}
    vwts::Vector{Int}
    hwts::Vector{Int}
end

struct __pindex__
    p1::Vector{Int}
    p2::Vector{Int}
end

struct __least_common_ancestor__
    rmq_sparse_table::Matrix{Int}
    euler_level::Vector{Int}
    child::Vector{Int}
    parents::Vector{Int}
    euler_tour::Vector{Int}
    level_vec::Vector{Int}
    fts::Vector{Int}
    ifts::Vector{Int}
end

struct __cut_profile__
    vtx_cuts::AbstractArray{Int}
    edge_cuts::Vector{Int}
    edge_diff::Vector{Int}
    pred::Vector{Int}
    edge_terminators::Vector{Int}
    p::__pindex__
    forced_type::Vector{Int}
    forced_0::Vector{Int}
    forced_1::Vector{Int}
    forced_01::Vector{Int}
    FB0::Vector{Int}
    FB1::Vector{Int}
    edge_cuts_0::Vector{Int}
    edge_cuts_1::Vector{Int}
end

struct __best_partition__
    total_cost::Float64
    area_cost::Float64
    cut_cost::Float64
    cutsize::Int
    partition::Vector{Int}
    cut_point::Int
    area::Vector{Int}
end

struct __recursive_parts__
    hypergraph::__hypergraph__
    T::SimpleWeightedGraphs.SimpleGraph
    distilled_cuts::__cut_profile__
    capacities::Vector{Int}
    cluster_labels::Vector{Int}
end

mutable struct __tree_cuts__
    nforced0::Int
    nforced1::Int
    nforced01::Int
    total_vwts::Int
    exc0::Vector{Int}
    exc1::Vector{Int}
    area::Vector{Int}
    cut_cost0::Vector{Float64}
    cut_cost1::Vector{Float64}
    ratio_cost::Vector{Float64}
    cut_cost::Vector{Float64}
    area_cost::Vector{Float64}
    total_cost::Vector{Float64}
    polarity::Vector{Int}
    status_flag::Vector{Int}
    area_util0::Vector{Int}
    area_util1::Vector{Int}
    pred::Vector{Int}
    hyperedges_flag::Vector{Int}
end

# Set KSPECPART_SERIAL_KERNELS=1 to force the serial code paths in the
# matrix-free operators. This makes runs bit-for-bit reproducible regardless of
# the thread count (parallel float reduction is deterministic for a fixed thread
# count but not identical to the serial summation order).
const SERIAL_KERNELS = get(ENV, "KSPECPART_SERIAL_KERNELS", "0") == "1"

# Below this hyperedge count `hypl` runs serially. `hypl` is memory-bandwidth
# bound and per-call cheap, so the task-spawn + per-thread-buffer overhead only
# pays off for large hypergraphs; benchmarking showed mid-size instances (~10^5
# hyperedges) are neutral-to-slower, so the threshold is set conservatively.
const HYPL_PARALLEL_THRESHOLD = 200_000

# Partition 1:n into (at most) k contiguous ranges. Fixed assignment so the
# parallel reduction order is deterministic.
function chunk_ranges(n::Int, k::Int)
    k = max(1, min(k, n))
    base, extra = divrem(n, k)
    ranges = Vector{UnitRange{Int}}(undef, k)
    start = 1
    @inbounds for t in 1:k
        len = base + (t <= extra ? 1 : 0)
        ranges[t] = start:(start + len - 1)
        start += len
    end
    return ranges
end

# Accumulate the hyperedge-Laplacian contribution of hyperedges `jrange` into y.
# The arithmetic is identical to the original serial loop, so the serial path
# (jrange == 1:m) is bit-for-bit unchanged.
function hypl_kernel!(y::Vector{Float64}, hypergraph::__hypergraph__,
                      x::AbstractArray, epsilon::Int, jrange::UnitRange{Int})
    eind = hypergraph.eind
    eptr = hypergraph.eptr
    w = hypergraph.hwts
    @inbounds for j in jrange
        first_valid_entry = eptr[j]
        first_invalid_entry = eptr[j+1]
        k = first_invalid_entry - first_valid_entry
        scale = (floor(k/2) * ceil(k/2))/(k-1)
        sm = 0.0
        for t in first_valid_entry:first_invalid_entry-1
            sm += x[eind[t]]
        end
        sm /= k
        for t in first_valid_entry:first_invalid_entry-1
            idx = eind[t]
            y[idx] += w[j] * (x[idx] - sm)/(scale*epsilon)
        end
    end
    return y
end

function hypl(hypergraph::__hypergraph__, x::AbstractArray, epsilon::Int)
    n = length(x)
    m = hypergraph.num_hyperedges
    nt = Threads.nthreads()
    if SERIAL_KERNELS || nt == 1 || m < HYPL_PARALLEL_THRESHOLD
        return hypl_kernel!(zeros(Float64, n), hypergraph, x, epsilon, 1:m)
    end
    # Each task accumulates into its own buffer (no write conflicts); buffers are
    # then reduced in a fixed order so the result is deterministic.
    ranges = chunk_ranges(m, nt)
    nb = length(ranges)
    ybuf = [zeros(Float64, n) for _ in 1:nb]
    @sync for t in 1:nb
        Threads.@spawn hypl_kernel!(ybuf[t], hypergraph, x, epsilon, ranges[t])
    end
    y = ybuf[1]
    @inbounds for t in 2:nb
        y .+= ybuf[t]
    end
    return y
end

function clique(x::AbstractArray, 
                vwts::Vector{Int}, 
                multiplier::AbstractArray)
    twt = sum(vwts)
    n = size(x, 1)
    y = zeros(Float64, n)
    s = multiplier[1]/twt
    kvec = vwts'x
    # O(n) and applied inside the (possibly block-parallel) eigensolves; kept
    # serial so it composes without nested-@threads oversubscription. Each y[j]
    # is independent, so this is numerically identical to the threaded version.
    @inbounds for j in 1:n
        y[j] += twt * ((vwts[j] * x[j]) - ((kvec * vwts[j])/twt)) * s
    end
    return y
end

function bi_clique(x::AbstractArray, pindex::__pindex__)
    n = length(x)
    y = zeros(n)
    d1 = ones(length(pindex.p1))
    d2 = ones(length(pindex.p2))
    t1 = Threads.@spawn (sum(d2) .* d1 .* x[pindex.p1] - d1 * (d2' * x[pindex.p2]))
    t2 = Threads.@spawn (sum(d1) .* d2 .* x[pindex.p2] - d2 * (d1' * x[pindex.p1])) 
    t1 = fetch(t1)
    t2 = fetch(t2)
    y[pindex.p1] = t1
    y[pindex.p2] = t2
    return y
end

function process_hint(partition::Vector{Int}, 
                    new_indices::Vector{Int},
                    processed_partition::Vector{Int})
    for i in 1:length(new_indices)
        mapped_i = new_indices[i]
        if mapped_i == 0
            continue
        end
        processed_partition[mapped_i] = partition[i]
    end
end 

function find_labels(clusters::AbstractArray, 
                    n::Int)
    labels = zeros(Int, n)
    for i in 1:length(clusters) 
        for j in 1:length(clusters[i])
            labels[clusters[i][j]] = i-1
        end
    end
    return labels
end

function write_partition(partition::Vector{Int},
                        partition_file_name::String)
    f = open(partition_file_name, "w")
    for i in 1:length(partition) 
        println(f, partition[i])
    end
    close(f)
end

# ---------------------------------------------------------------------------
# External-tool and scratch-directory configuration.
#
# Defaults are derived from the directory containing this file, so K-SpecPart
# runs from wherever it is checked out. Every path can be overridden with an
# environment variable -- the recommended way to point at your own binaries
# without editing source.
#
#   KSPECPART_SOURCE_DIR  scratch directory for temporary files (default: here)
#   KSPECPART_METIS_DIR   directory holding metis_script.sh   (default: here)
#   KSPECPART_HMETIS      hMETIS binary
#   KSPECPART_ILP         OR-Tools / CPLEX ILP partitioner binary
#   KSPECPART_REFINER     TritonPart / OpenROAD FM refiner binary
# ---------------------------------------------------------------------------

const _PKG_DIR = @__DIR__

# A trailing space is intentionally appended to the binary paths because every
# caller builds command strings by concatenating arguments after the path.
source_dir = get(ENV, "KSPECPART_SOURCE_DIR", _PKG_DIR)
metis_path = get(ENV, "KSPECPART_METIS_DIR", _PKG_DIR)
hmetis_path = get(ENV, "KSPECPART_HMETIS",
                  "/home/fetzfs_projects/SpecPart/K_SpecPart/hmetis") * " "
# Optional direct k-way partitioner (khmetis / HMETIS_PartKway). When set, it is
# used for the overlay solve at k > 2 with a direct per-block imbalance bound
# (UBfactor = k*epsilon), as recommended by the hMETIS manual for large k. Empty
# => fall back to recursive-bisection hmetis. Set via KSPECPART_KHMETIS.
khmetis_path = get(ENV, "KSPECPART_KHMETIS", "")
ilp_path = get(ENV, "KSPECPART_ILP",
               joinpath(_PKG_DIR, "ilp_partitioner", "build", "ilp_part")) * " "
triton_part_refiner_path = get(ENV, "KSPECPART_REFINER",
                               "/home/bodhi91/TritonPart_OpenROAD/build/src/openroad") * " "

# The OpenROAD refiner and the ILP partitioner are dynamically linked against
# OR-Tools (libortools.so.9), which is found only via LD_LIBRARY_PATH (no rpath).
# Prepend its directory to LD_LIBRARY_PATH here so every child process inherits
# it regardless of the shell the run was launched from. Override the directory
# with KSPECPART_ORTOOLS_LIB (set it empty to disable this behavior).
const ORTOOLS_LIB = get(ENV, "KSPECPART_ORTOOLS_LIB",
    "/home/tool/ortools/install/or-tools_cpp_CentOSStream-8-64bit_v9.4.1874/lib64")
if !isempty(ORTOOLS_LIB)
    existing = get(ENV, "LD_LIBRARY_PATH", "")
    if !occursin(ORTOOLS_LIB, existing)
        ENV["LD_LIBRARY_PATH"] = isempty(existing) ? ORTOOLS_LIB : ORTOOLS_LIB * ":" * existing
    end
end
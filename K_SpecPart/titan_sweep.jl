# Titan23 multi-seed sweep for K-SpecPart.
#
# Runs every Titan design x K in {2,3,4} x several seeds, starting from the
# provided hint solutions, and records for each run: the hint cut, the output
# cut, per-block balance, and wall-clock time. Every output partition is saved
# to disk.
#
# Usage (from the K_SpecPart directory):
#     julia -t 8 --project=. titan_sweep.jl
#
# - Results CSV is appended after every run, so partial progress is never lost.
# - The sweep is RESUMABLE: a (design, K, seed) whose partition file already
#   exists is skipped. Delete a partition file (or the whole output dir) to
#   force a rerun.
# - Each run is wrapped in try/catch; one failure does not abort the sweep.
#
# Scope can be narrowed via environment variables (handy for testing / sharding):
#     SWEEP_DESIGNS=neuron,dart  SWEEP_KS=2,4  SWEEP_SEEDS=0  \
#         SWEEP_OUT=/tmp/sweep  julia -t 8 --project=. titan_sweep.jl
#
# Tip: redirect the verbose solver logs, e.g.
#     julia -t 8 --project=. titan_sweep.jl > sweep.stdout.log 2>&1

include("specpart.jl")
using .SpecPart
include("collect_best_cuts.jl")

# --------------------------------------------------------------------------
# Configuration (edit as needed)
# --------------------------------------------------------------------------
const BENCH_DIR = "/home/fetzfs_projects/SpecPart/HypergraphPartitioning/benchmark/Titan23_benchmark"
const HINT_DIR  = "/home/fetzfs_projects/SpecPart/HypergraphPartitioning/K_specpart_solutions/Titan23_benchmarks/ub_factor_2"
const OUT_DIR   = get(ENV, "SWEEP_OUT", joinpath(@__DIR__, "titan_sweep_results"))

parse_int_list(s, default) = isempty(strip(s)) ? default :
    [parse(Int, strip(x)) for x in split(s, ",")]

const KS    = parse_int_list(get(ENV, "SWEEP_KS", ""), [2, 3, 4])      # K values
const SEEDS = parse_int_list(get(ENV, "SWEEP_SEEDS", ""), [0, 1, 2])   # seeds per (design,K)
const IMB   = 2                    # imbalance factor (matches ub_factor_2 hints)

# K-SpecPart parameters (these mirror the current shipped defaults).
const EIGVECS      = 2
const REFINE_ITERS = 2
const SOLVER_ITERS = 40
const BEST_SOLNS   = 5
const NCYCLES      = 1

# Optionally restrict the design list; empty => all *.hgr in BENCH_DIR.
# Override with SWEEP_DESIGNS=design1,design2,...
const ONLY_DESIGNS = isempty(strip(get(ENV, "SWEEP_DESIGNS", ""))) ? String[] :
    [strip(x) for x in split(ENV["SWEEP_DESIGNS"], ",")]
# --------------------------------------------------------------------------

hint_path(design, k) =
    joinpath(HINT_DIR, string(k) * "_way",
             design * ".hgr.specpart.ubfactor." * string(IMB) * ".part." * string(k))

part_out_path(design, k, seed) =
    joinpath(OUT_DIR, "partitions", "$(design).k$(k).seed$(seed).part")

function discover_designs()
    if !isempty(ONLY_DESIGNS)
        return sort(ONLY_DESIGNS)
    end
    files = filter(f -> endswith(f, ".hgr"), readdir(BENCH_DIR))
    return sort([replace(f, r"\.hgr$" => "") for f in files])
end

function append_csv_row(csv_file, row)
    open(csv_file, "a") do io
        println(io, join(row, ","))
    end
end

function main()
    SpecPart.setup_logging()
    mkpath(OUT_DIR)
    mkpath(joinpath(OUT_DIR, "partitions"))
    csv_file = joinpath(OUT_DIR, "sweep_results.csv")
    if !isfile(csv_file)
        append_csv_row(csv_file,
            ["design", "num_parts", "seed", "hint_cut", "out_cut",
             "balance", "runtime_s", "status", "partition_file"])
    end

    designs = discover_designs()
    @info "Sweep: $(length(designs)) designs x $(length(KS)) K x $(length(SEEDS)) seeds = " *
          "$(length(designs)*length(KS)*length(SEEDS)) runs; threads=$(Threads.nthreads())"

    for design in designs
        hgr = joinpath(BENCH_DIR, design * ".hgr")
        if !isfile(hgr)
            @warn "missing hypergraph, skipping" design
            continue
        end
        # Read the hypergraph once per design for hint-cut / balance reporting.
        hg = SpecPart.read_hypergraph_file(hgr)
        for k in KS
            hint = hint_path(design, k)
            if !isfile(hint)
                @warn "missing hint, skipping" design k hint
                continue
            end
            # Hint cut for reference (full-vertex indexing).
            hint_cut = -1
            try
                hp = SpecPart.read_hint_file(hint)
                (hint_cut, _) = SpecPart.golden_evaluator(hg, k, hp)
            catch e
                @warn "could not evaluate hint cut" design k exception=e
            end
            for seed in SEEDS
                pfile = part_out_path(design, k, seed)
                if isfile(pfile)
                    @info "skip (exists)" design k seed
                    continue
                end
                @info "[sweep] RUN design=$design K=$k seed=$seed hint_cut=$hint_cut"
                status = "OK"
                out_cut = -1
                balance = ""
                runtime = -1.0
                part = Int[]
                try
                    t = @elapsed begin
                        (part, out_cut) = SpecPart.specpart_run(hgr;
                            hint_file   = hint,
                            imb         = IMB,
                            num_parts   = k,
                            eigvecs     = EIGVECS,
                            refine_iters= REFINE_ITERS,
                            solver_iters= SOLVER_ITERS,
                            best_solns  = BEST_SOLNS,
                            ncycles     = NCYCLES,
                            seed        = seed)
                    end
                    runtime = round(t, digits=2)
                    # Balance is only meaningful if the returned partition has
                    # full length (no islands were removed); otherwise leave blank.
                    if length(part) == hg.num_vertices
                        (_, bal) = SpecPart.golden_evaluator(hg, k, part)
                        balance = join(bal, "|")
                    end
                    SpecPart.write_partition(part, pfile)
                catch e
                    status = "ERROR"
                    @error "run failed" design k seed exception=(e, catch_backtrace())
                end
                append_csv_row(csv_file,
                    [design, k, seed, hint_cut, out_cut, balance, runtime,
                     status, isfile(pfile) ? pfile : ""])
                @info "[sweep] DONE design=$design K=$k seed=$seed out_cut=$out_cut " *
                      "runtime_s=$runtime status=$status"
            end
        end
    end
    @info "Sweep complete. Results: $csv_file ; partitions in $(joinpath(OUT_DIR, "partitions"))"
    # Aggregate the best cut per (design, K) into summary CSVs.
    try
        collect_best_cuts(OUT_DIR)
    catch e
        @warn "best-cut collection failed" exception=e
    end
end

main()

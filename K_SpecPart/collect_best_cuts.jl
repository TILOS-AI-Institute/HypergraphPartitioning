# Collect the best cuts from a titan_sweep run into summary CSVs.
#
# Reads  <OUT_DIR>/sweep_results.csv  (written by titan_sweep.jl) and produces:
#   <OUT_DIR>/best_cuts.csv        long form, one row per (design, K):
#                                  design, num_parts, hint_cut, best_cut,
#                                  best_seed, improvement_pct, mean_cut,
#                                  num_seeds, mean_runtime_s, partition_file
#   <OUT_DIR>/best_cuts_pivot.csv  wide form: design, best_k2, best_k3, best_k4, ...
#
# Usage:
#   julia collect_best_cuts.jl            # uses ./titan_sweep_results
#   SWEEP_OUT=/path julia collect_best_cuts.jl
# It can be run any time, including while a sweep is still in progress.

using Printf
using Statistics

# One aggregated record per (design, K).
struct BestRecord
    design::String
    k::Int
    hint_cut::Int
    best_cut::Int
    best_seed::Int
    mean_cut::Float64
    num_seeds::Int
    mean_runtime::Float64
    partition_file::String
end

function collect_best_cuts(out_dir::AbstractString)
    csv_in = joinpath(out_dir, "sweep_results.csv")
    isfile(csv_in) || error("results file not found: $csv_in")

    # group[(design, k)] = Vector of (seed, hint_cut, out_cut, runtime, pfile)
    groups = Dict{Tuple{String,Int}, Vector{NTuple{5,Any}}}()
    open(csv_in, "r") do io
        for (lineno, ln) in enumerate(eachline(io))
            lineno == 1 && continue                  # header
            isempty(strip(ln)) && continue
            f = split(ln, ',')
            length(f) < 8 && continue
            status = strip(f[8])
            status == "OK" || continue
            out_cut = tryparse(Int, strip(f[5]))
            (out_cut === nothing || out_cut < 0) && continue
            design = String(strip(f[1]))
            k = parse(Int, strip(f[2]))
            seed = parse(Int, strip(f[3]))
            hint_cut = something(tryparse(Int, strip(f[4])), -1)
            runtime = something(tryparse(Float64, strip(f[7])), -1.0)
            pfile = length(f) >= 9 ? String(strip(f[9])) : ""
            push!(get!(groups, (design, k), Vector{NTuple{5,Any}}()),
                  (seed, hint_cut, out_cut, runtime, pfile))
        end
    end

    records = BestRecord[]
    for ((design, k), rows) in groups
        cuts = [r[3] for r in rows]
        best_i = argmin(cuts)
        hint = maximum(r[2] for r in rows)           # constant across seeds; -1 if unknown
        push!(records, BestRecord(design, k, hint, cuts[best_i], rows[best_i][1],
            mean(cuts), length(rows), mean(r[4] for r in rows), rows[best_i][5]))
    end
    sort!(records, by = r -> (r.design, r.k))

    # Long-form summary.
    best_csv = joinpath(out_dir, "best_cuts.csv")
    open(best_csv, "w") do io
        println(io, "design,num_parts,hint_cut,best_cut,best_seed,improvement_pct," *
                    "mean_cut,num_seeds,mean_runtime_s,partition_file")
        for r in records
            imp = r.hint_cut > 0 ?
                @sprintf("%.2f", 100 * (r.hint_cut - r.best_cut) / r.hint_cut) : "NA"
            println(io, join((r.design, r.k, r.hint_cut, r.best_cut, r.best_seed,
                imp, @sprintf("%.2f", r.mean_cut), r.num_seeds,
                @sprintf("%.2f", r.mean_runtime), r.partition_file), ","))
        end
    end

    # Wide pivot: design x K best cut.
    ks = sort(unique(r.k for r in records))
    by_design = Dict{String, Dict{Int,Int}}()
    for r in records
        get!(by_design, r.design, Dict{Int,Int}())[r.k] = r.best_cut
    end
    pivot_csv = joinpath(out_dir, "best_cuts_pivot.csv")
    open(pivot_csv, "w") do io
        println(io, "design," * join(["best_k$(k)" for k in ks], ","))
        for design in sort(collect(keys(by_design)))
            vals = [haskey(by_design[design], k) ? string(by_design[design][k]) : ""
                    for k in ks]
            println(io, design * "," * join(vals, ","))
        end
    end

    @info "Wrote $(length(records)) (design,K) best-cut records"
    @info "  long : $best_csv"
    @info "  pivot: $pivot_csv"
    return best_csv, pivot_csv
end

if abspath(PROGRAM_FILE) == @__FILE__
    out_dir = get(ENV, "SWEEP_OUT", joinpath(@__DIR__, "titan_sweep_results"))
    collect_best_cuts(out_dir)
end

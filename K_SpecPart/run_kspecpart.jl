#!/usr/bin/env julia
#
# Command-line front-end for K-SpecPart.
#
# Usage:
#   julia --project=. -t auto run_kspecpart.jl <hypergraph.hgr> [options]
#
# Required:
#   <hypergraph.hgr>        input hypergraph in hMETIS format
#
# Options:
#   --hint FILE             initial hint partition (hMETIS .part.K format).
#                           Strongly recommended: K-SpecPart *refines* a hint.
#   --fixed FILE            fixed-vertex file
#   --k K                   number of blocks           (default 2)
#   --imb PCT               per-block imbalance in %    (default 2)
#   --eigvecs N             eigenvectors per problem    (default 2)
#   --refine-iters N        refinement iterations       (default 2)
#   --solver-iters N        LOBPCG iterations           (default 40)
#   --best-solns N          candidates kept for overlay (default 5)
#   --ncycles N             sparsifier cycles           (default 1)
#   --seed N                RNG seed                    (default 0)
#   --out FILE              also write the partition (one block id per line)
#   -h, --help              show this help and exit
#
# External tool paths are configured via environment variables (see README:
# "Configuration"): KSPECPART_HMETIS, KSPECPART_ILP, KSPECPART_REFINER,
# GPMETIS, KSPECPART_KHMETIS, KSPECPART_ORTOOLS_LIB, KSPECPART_SOURCE_DIR.

function usage()
    print("""
    K-SpecPart -- supervised spectral hypergraph partitioner

    Usage:
      julia --project=. -t auto run_kspecpart.jl <hypergraph.hgr> [options]

    Options:
      --hint FILE          initial hint partition (hMETIS .part.K); recommended
      --fixed FILE         fixed-vertex file
      --k K                number of blocks           (default 2)
      --imb PCT            per-block imbalance in %    (default 2)
      --eigvecs N          eigenvectors per problem    (default 2)
      --refine-iters N     refinement iterations       (default 2)
      --solver-iters N     LOBPCG iterations           (default 40)
      --best-solns N       candidates kept for overlay (default 5)
      --ncycles N          sparsifier cycles           (default 1)
      --seed N             RNG seed                    (default 0)
      --out FILE           write the partition (one block id per line)
      -V, --version        print version and exit
      -h, --help           show this help and exit
    """)
end

function version_string()
    try
        return strip(read(joinpath(@__DIR__, "VERSION"), String))
    catch
        return "unknown"
    end
end

function parse_args(args::Vector{String})
    if !isempty(args) && args[1] in ("--version", "-V")
        println("K-SpecPart ", version_string())
        exit(0)
    end
    if isempty(args) || args[1] in ("-h", "--help")
        usage()
        exit(isempty(args) ? 1 : 0)
    end
    hgr = args[1]
    opts = Dict{String,Any}(
        "hint" => "", "fixed" => "", "k" => 2, "imb" => 2.0, "eigvecs" => 2,
        "refine-iters" => 2, "solver-iters" => 40, "best-solns" => 5,
        "ncycles" => 1, "seed" => 0, "out" => "",
    )
    i = 2
    while i <= length(args)
        flag = args[i]
        startswith(flag, "--") || error("unexpected argument: $flag (use --help)")
        key = flag[3:end]
        haskey(opts, key) || error("unknown option: $flag (use --help)")
        i + 1 <= length(args) || error("missing value for $flag")
        val = args[i + 1]
        opts[key] = key in ("hint", "fixed", "out") ? val :
                    key == "imb" ? parse(Float64, val) : parse(Int, val)
        i += 2
    end
    return hgr, opts
end

const HGR, OPTS = parse_args(copy(ARGS))

include(joinpath(@__DIR__, "specpart.jl"))

partition, cutsize = Main.SpecPart.specpart_run(
    HGR;
    hint_file             = OPTS["hint"],
    hypergraph_fixed_file = OPTS["fixed"],
    num_parts             = OPTS["k"],
    imb                   = OPTS["imb"],
    eigvecs               = OPTS["eigvecs"],
    refine_iters          = OPTS["refine-iters"],
    solver_iters          = OPTS["solver-iters"],
    best_solns            = OPTS["best-solns"],
    ncycles               = OPTS["ncycles"],
    seed                  = OPTS["seed"],
)

@info "K-SpecPart finished" cutsize=cutsize k=OPTS["k"]

if !isempty(OPTS["out"])
    open(OPTS["out"], "w") do io
        for p in partition
            println(io, p)
        end
    end
    @info "Partition written" file=OPTS["out"]
end

function write_hypergraph(hgraph::__hypergraph__, fname::String)
    n = hgraph.num_vertices
    e = hgraph.num_hyperedges
    hedges = hgraph.eind
    eptr = hgraph.eptr
    hwts = hgraph.hwts
    vwts = hgraph.vwts
    fixed = hgraph.fixed
    wt_flag = maximum(vwts) > 1 ? true : false
    f = open(fname, "w")
    println(f, e, " ", n, " 11")

    for i in 1:e
        start_idx = eptr[i]
        end_idx = eptr[i+1]-1
        print(f, hwts[i])
        for j in start_idx:end_idx
            print(f, " ", hedges[j])
        end
        print(f, "\n")
    end
    for i in 1:n
        if wt_flag == 0
            println(f, vwts[i])
        else
            println(f, vwts[i]+1)
        end
    end

    close(f)

    if maximum(fixed) > -1
        f = open(fname*".fixed", "w")
        for i in 1:n
            # Was `fixed_vtxs[i]`, an undefined variable that would throw
            # whenever the hypergraph actually had fixed vertices.
            println(f, fixed[i])
        end
        close(f)
    end
end

function check_balance(hgraph::__hypergraph__, partition::Vector{Int}, num_parts::Int, ub_factor::Real)
    blocks = zeros(Int, num_parts)
    for i in 1:length(partition)
        blocks[partition[i]+1] += hgraph.vwts[i]
    end

    # Per-block capacity matching the internal capacities: (100/k + ub_factor)%.
    # (Previously used the 2-way (50 + ub_factor)% formula regardless of k.)
    max_balance = Int(ceil((100/num_parts + ub_factor) * sum(hgraph.vwts)/100))
    for i in 1:num_parts
        if blocks[i] > max_balance
            return false
        end
    end
    return true
end

# Read a partition file (one block id per line) into `partition`, in place.
function read_partition_file!(partition::Vector{Int}, pfile::String)
    itr = 0
    open(pfile, "r") do f
        for ln in eachline(f)
            itr += 1
            partition[itr] = parse(Int, ln)
        end
    end
    return partition
end

use_khmetis(num_parts::Int) = num_parts > 2 && !isempty(strip(khmetis_path))

# Convert K-SpecPart's per-block percentage-point slack `ub_factor` (each block
# is allowed up to (100/k + ub_factor)% of the total weight) into a VALID
# INTEGER hMETIS UBfactor. The mapping depends on the program (see hMETIS
# manual):
#   - k = 2 (recursive bisection): UBfactor = ub_factor   -> (50 +- ub_factor)%.
#   - k > 2 via khmetis (direct k-way): UBfactor = k*ub_factor (the paper's
#     "k*epsilon"); heaviest block <= (1 + UBfactor/100)*(W/k).
#   - k > 2 via recursive hmetis: UBfactor is applied PER BISECTION and
#     compounds over ceil(log2 k) levels, so we invert that to approximate the
#     target per-block tolerance. (For large/odd k, khmetis is recommended.)
# The result is clamped to hMETIS's valid range [1, 49].
function balance_ubfactor(num_parts::Int, ub_factor::Real)
    eps_rel = num_parts * ub_factor / 100        # heaviest block <= (1+eps_rel)*W/k
    if num_parts == 2
        b = round(Int, ub_factor)
    elseif use_khmetis(num_parts)
        b = round(Int, 100 * eps_rel)            # == round(k * ub_factor)
    else
        levels = max(1, ceil(Int, log2(num_parts)))
        b = round(Int, 50 * ((1 + eps_rel)^(1/levels) - 1))
    end
    return clamp(b, 1, 49)
end

# Full shell command for the overlay-solve partitioner: khmetis (direct k-way)
# when available and k>2, otherwise the recursive-bisection hmetis.
function partitioner_command_string(hmetis_path::String, hgr_file::String,
                                    num_parts::Int, ub_factor::Real)
    ub = balance_ubfactor(num_parts, ub_factor)
    runs = 10; ctype = 1; vcycle = 1; dbglvl = 0
    if use_khmetis(num_parts)
        otype = 1   # minimize hyperedge cut
        return khmetis_path * " " * hgr_file * " " * string(num_parts) * " " *
               string(ub) * " " * string(runs) * " " * string(ctype) * " " *
               string(otype) * " " * string(vcycle) * " " * string(dbglvl)
    else
        rtype = 1; reconst = 0
        return hmetis_path * " " * hgr_file * " " * string(num_parts) * " " *
               string(ub) * " " * string(runs) * " " * string(ctype) * " " *
               string(rtype) * " " * string(vcycle) * " " * string(reconst) * " " *
               string(dbglvl)
    end
end

# Exact-ILP size limits for the overlay-contracted hypergraph. Raised from the
# original (1500 / 300) because the contracted hypergraph is small by
# construction; the wall-clock limit below bounds the worst case so larger
# instances that the ILP cannot crack in time fall back to hMETIS.
# Set KSPECPART_ILP_BOOST=0 to recover the original thresholds.
const ILP_BOOST = get(ENV, "KSPECPART_ILP_BOOST", "1") == "1"
const ILP_MAX_HE_2WAY = ILP_BOOST ? 3000 : 1500
const ILP_MAX_HE_KWAY = ILP_BOOST ? 800 : 300
const ILP_TIME_LIMIT_S = 20

# Run an external command with a wall-clock limit; kill it on timeout. Returns
# true only if it finished on its own and exited successfully.
function run_with_timeout(cmd::Cmd, seconds::Real)
    p = run(cmd, wait=false)
    timedout = Ref(false)
    timer = Timer(seconds) do _
        if process_running(p)
            timedout[] = true
            kill(p, 9)
        end
    end
    wait(p)
    close(timer)
    return !timedout[] && success(p)
end

function optimal_partitioner(hmetis_path::String, cplex_path::String, hgraph::__hypergraph__, num_parts::Int, ub_factor::Real)
    partition = zeros(Int, hgraph.num_vertices)
    hgr_file_name = source_dir * "/" * "coarse.hgr"
    write_hypergraph(hgraph, hgr_file_name)
    # The hypergraph is small enough for the exact ILP solver (the two former
    # branches had identical bodies and only differed in this size guard).
    use_ilp = (num_parts == 2 && hgraph.num_hyperedges < ILP_MAX_HE_2WAY) ||
              (num_parts > 2 && hgraph.num_hyperedges < ILP_MAX_HE_KWAY)
    if use_ilp
        ilp_string = ilp_path * " " * hgr_file_name * " " * string(num_parts) * " " * string(ub_factor)
        pfile = hgr_file_name * ".part." * string(num_parts)
        # Time-bounded exact solve; if it cannot finish in time (or returns an
        # infeasible/degenerate result), fall back to hMETIS.
        ilp_done = run_with_timeout(`sh -c $ilp_string`, ILP_TIME_LIMIT_S)
        ilp_ok = false
        if ilp_done && isfile(pfile)
            read_partition_file!(partition, pfile)
            (cutsize, ~) = golden_evaluator(hgraph, num_parts, partition)
            ilp_ok = check_balance(hgraph, partition, num_parts, ub_factor) && cutsize != 0
        end
        if !ilp_ok
            hmetis_string = partitioner_command_string(hmetis_path, hgr_file_name, num_parts, ub_factor)
            run(`sh -c $hmetis_string`, wait=true)
            read_partition_file!(partition, pfile)
        end
        run(`rm -f $hgr_file_name`, wait=true)
        run(`rm -f $pfile`, wait=true)
    else
        # 10 parallel runs of hMETIS; keep the lowest-cut result.
        parallel_runs = 10
        partitions = [zeros(Int, length(partition)) for i in 1:parallel_runs]
        cutsizes = zeros(Int, parallel_runs)
        @sync Threads.@threads for i in 1:parallel_runs
            local_hgr_name = hgr_file_name * "." * string(i)
            run(`sh -c $("cp " * hgr_file_name * " " * local_hgr_name)`, wait=true)
            hmetis_string = partitioner_command_string(hmetis_path, local_hgr_name, num_parts, ub_factor)
            run(`sh -c $hmetis_string`, wait=true)
        end
        for i in 1:parallel_runs
            local_hgr_name = hgr_file_name * "." * string(i)
            local_pfile_name = local_hgr_name * ".part." * string(num_parts)
            read_partition_file!(partitions[i], local_pfile_name)
            (cutsizes[i], ~) = golden_evaluator(hgraph, num_parts, partitions[i])
            run(`rm $local_hgr_name $local_pfile_name`, wait=true)
        end
        ~, best_cut_idx = findmin(cutsizes)
        run(`rm $hgr_file_name`, wait=true)
        partition = partitions[best_cut_idx]
    end
    return partition
end
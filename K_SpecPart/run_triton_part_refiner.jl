function triton_part_refine(refiner_path::String, hypergraph::String, partition::String, num_parts::Int, ub_factor::Real, seed::Int, id::Int)
    # `triton_part_refine` maps `-balance_constraint` directly onto its UBfactor
    # (TritonPart.cpp: ub_factor_ = balance_constraint_arg), and it defaults to
    # 1.0. Passing the actual imbalance gives FM the legal slack it needs for
    # imb > 1 instead of being over-constrained. `-ub_factor` is NOT a valid
    # flag (parse_key_args errors on unknown keys), so the imbalance must go
    # through `-balance_constraint`. Set KSPECPART_FM_BALANCE=0 to omit it
    # (recovering the old default UBfactor of 1.0).
    balance_flag = get(ENV, "KSPECPART_FM_BALANCE", "1") == "1" ?
                   " -balance_constraint " * string(ub_factor) : ""
    line = "triton_part_refine" *
            " -hypergraph_file " * hypergraph *
            " -partition_file " * partition *
            " -num_parts " * string(num_parts) *
            balance_flag *
            " -seed " * string(seed)
    tcl_name = joinpath(source_dir, "run_triton_part_refiner." * string(id) * ".tcl")
    log_name = "run.log." * string(id)
    sh_name = joinpath(source_dir, "run_refiner." * string(id) * ".sh")
    f = open(tcl_name, "w")
    #f = open(source_dir * "run_triton_part_refiner.tcl", "w")
    println(f, line) 
    println(f, "exit")
    close(f)
    f = open(sh_name, "w")
    # Use bash with pipefail so a refiner failure (e.g. a missing shared library)
    # propagates as a non-zero exit instead of being masked by `tee`'s success.
    # Previously the masked failure silently skipped FM refinement and degraded
    # results without any error.
    println(f, "#!/bin/bash")
    println(f, "set -o pipefail")
    println(f, refiner_path * tcl_name * " | tee " * log_name)
    close(f)
    cmd = "chmod 777 " * sh_name 
    run(`sh -c $cmd`, wait=true)
    cmd = "chmod 777 " * tcl_name
    run(`sh -c $cmd`, wait = true)
    triton_part_refiner_cmd = `$sh_name
                                $hypergraph
                                $partition
                                $num_parts
                                $ub_factor
                                $seed`
    #=triton_part_refiner_cmd = `/home/fetzfs_projects/SpecPart/src/run_fm_refinement.sh
                                    $hypergraph
                                    $partition
                                    $num_parts
                                    $ub_factor
                                    $seed`=#
    run(triton_part_refiner_cmd, wait=true) 
    rm = "rm " * tcl_name * " " * log_name * " " * sh_name
    run(`sh -c $rm`, wait = true)
end
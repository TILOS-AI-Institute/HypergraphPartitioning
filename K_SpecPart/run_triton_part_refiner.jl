function triton_part_refine(refiner_path::String, hypergraph::String, partition::String, num_parts::Int, ub_factor::Int, seed::Int, id::Int)
    line = "triton_part_refine" * 
            " -hypergraph_file " * hypergraph *
            " -partition_file " * partition *
            " -num_parts " * string(num_parts)
            " -ub_factor " * string(ub_factor)
            " -seed " * string(seed)
    tcl_name = source_dir * "run_triton_part_refiner." * string(id) * ".tcl"
    log_name = "run.log." * string(id)
    sh_name = source_dir * "run_refiner." * string(id) * ".sh"
    f = open(tcl_name, "w")
    #f = open(source_dir * "run_triton_part_refiner.tcl", "w")
    println(f, line) 
    println(f, "exit")
    close(f)
    f = open(sh_name, "w")
    line = refiner_path * 
                tcl_name * " | tee " * log_name
    println(f, line)
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
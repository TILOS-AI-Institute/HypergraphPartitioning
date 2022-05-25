include("SpectralRefinement.jl")    

using Random
using Statistics

Random.seed!(1)

csv_file = "default_hmetis_rpt_julia_titan_500_runs/hmetis_Titan_summary_file.txt"

nsamples = 499
multi_starts = Vector{Int}(1:20)
samples = 100

gsm_seeds = Vector{Int}(1:nsamples)
gsm_cuts = zeros(Int, nsamples)
gsm_runtime = zeros(nsamples)

sparc_seeds = Vector{Int}(1:nsamples)
sparc_cuts = zeros(Int, nsamples)
sparc_runtime = zeros(nsamples)

f = open(csv_file)
global i = 0

for ln in eachline(f)
    r = match(r"([a-zA-Z_0-9]+),([0-9]+),([0-9]+),([0-9]+),([0-9]+),([0-9.]+)", ln)

    if r == nothing
        continue
    end

    r_captures = r.captures
    benchmark = r_captures[1]
    seed = parse(Int, r_captures[4])
    cut = parse(Int, r_captures[5])
    runtime = parse(Float64, r_captures[6])
    
    
    if benchmark == "gsm_switch"
        global i += 1
        gsm_cuts[i] = cut
        gsm_runtime[i] = runtime
        gsm_seeds[i] = seed
    else
        if i == nsamples
            global i = 0
        end

        global i += 1
        sparc_cuts[i] = cut
        sparc_runtime[i] = runtime
        sparc_seeds[i] = seed
    end
end

close(f)

opfile = "qor_vs_runtime.csv"

f = open(opfile, "w")
println(f, "Benchmark,Starts,Avg_hmetis_cut,Avg_hmetis_time,Avg_ilp_cut,Avg_ilp_time")
close(f)

for j in 1:length(multi_starts)
    best_hmetis_cut_list = zeros(Int, samples)
    best_hmetis_runtime = zeros(samples)
    ilp_cut_list = zeros(Int, samples)
    ilp_runtime = zeros(samples)
    multi_start = multi_starts[j]
    
    for k in 1:samples
        @info "Running multi-start #$multi_start with sample #$(k)"
        picked_samples = [rand(1:nsamples) for z in 1:multi_start]   
        picked_cuts = gsm_cuts[picked_samples]
        picked_runtime = gsm_runtime[picked_samples]
        picked_seeds = gsm_seeds[picked_samples]
        best_picked_cut, idx = findmin(picked_cuts)
        best_picked_runtime = picked_runtime[idx]
        best_hmetis_cut_list[k] += best_picked_cut
        best_hmetis_runtime[k] += best_picked_runtime

        if multi_start >= 5
            sperm = sortperm(picked_cuts)
            best_seeds = picked_seeds[sperm[1:5]]
            t_overlay = @elapsed ilp_cut, ilp_run = Main.SpectralRefinement.OverlayBasedClusteringAndILP("/home/bodhi91/sandbox/Partitioners/Benchmarks/Titan/gsm_switch.hgr", 10, best_seeds)
            ilp_cut_list[k] = ilp_cut
            ilp_runtime[k] = t_overlay
        end
    end

    avg_best_hmetis_cut = mean(best_hmetis_cut_list)
    avg_best_hmetis_runtime = mean(best_hmetis_runtime) * multi_start
    avg_ilp_cut = mean(ilp_cut_list)
    avg_ilp_runtime = mean(ilp_runtime) + avg_best_hmetis_runtime

    f = open(opfile, "a")
    line = "gsm_switch"
    line *= "," * string(multi_start)
    line *= "," * string(avg_best_hmetis_cut)
    line *= "," * string(avg_best_hmetis_runtime)
    line *= "," * string(avg_ilp_cut)
    line *= "," * string(avg_ilp_runtime)
    line *= "\n"

    write(f, line)
    close(f)
end

for j in 1:length(multi_starts)
    best_hmetis_cut_list = zeros(Int, samples)
    best_hmetis_runtime = zeros(samples)
    ilp_cut_list = zeros(Int, samples)
    ilp_runtime = zeros(samples)
    multi_start = multi_starts[j]
    
    for k in 1:samples
        @info "Running multi-start #$multi_start with sample #$(k)"
        picked_samples = [rand(1:nsamples) for z in 1:multi_start]   
        picked_cuts = sparc_cuts[picked_samples]
        picked_runtime = sparc_runtime[picked_samples]
        picked_seeds = sparc_seeds[picked_samples]
        best_picked_cut, idx = findmin(picked_cuts)
        best_picked_runtime = picked_runtime[idx]
        best_hmetis_cut_list[k] += best_picked_cut
        best_hmetis_runtime[k] += best_picked_runtime

        if multi_start >= 5
            sperm = sortperm(picked_cuts)
            best_seeds = picked_seeds[sperm[1:5]]
            ilp_cut, ilp_run = Main.SpectralRefinement.OverlayBasedClusteringAndILP("/home/bodhi91/sandbox/Partitioners/Benchmarks/Titan/sparcT2_core.hgr", 10, best_seeds)
            ilp_cut_list[k] += ilp_cut
            ilp_runtime[k] += ilp_run
        end
    end

    avg_best_hmetis_cut = mean(best_hmetis_cut_list)
    avg_best_hmetis_runtime = mean(best_hmetis_runtime) * multi_start
    avg_ilp_cut = mean(ilp_cut_list)
    avg_ilp_runtime = mean(ilp_runtime) + avg_best_hmetis_runtime

    f = open(opfile, "a")
    line = "sparcT2_core"
    line *= "," * string(multi_start)
    line *= "," * string(avg_best_hmetis_cut)
    line *= "," * string(avg_best_hmetis_runtime)
    line *= "," * string(avg_ilp_cut)
    line *= "," * string(avg_ilp_runtime)
    line *= "\n"

    write(f, line)
    close(f)
end


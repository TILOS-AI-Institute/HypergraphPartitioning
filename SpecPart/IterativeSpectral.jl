include("SpectralRefinement.jl")

hypergraph_file = ARGS[1]
partition_file = ARGS[2]
num_parts = parse(Int, ARGS[3])
num_cycles = parse(Int, ARGS[4])
h_threshold = parse(Int, ARGS[5])
ub_factor = parse(Int, ARGS[6])
eig_vecs = parse(Int, ARGS[7])
solver_iters = parse(Int, ARGS[8])
spec_iters = parse(Int, ARGS[9])
p_solns = parse(Int, ARGS[10])
seed = parse(Int, ARGS[11])

spec_cut = Main.SpectralRefinement.SpectralHmetisRefinement(hg = hypergraph_file, pfile = partition_file, Nparts = num_parts, cycles = num_cycles, hyperedges_threshold = h_threshold, 
                                                    ub = ub_factor, nev = eig_vecs, seed = seed, solver_iters=solver_iters, refine_iters=spec_iters, best_solns=p_solns);
                                               
println("Supervised spectral cut: ", spec_cut)

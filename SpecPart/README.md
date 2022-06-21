# SpecPart #

SpecPart is the first supervised spectral hypergraph partitioning framework which enhances a given hypergraph partition solution.
This directory contains the Julia implementation of SpecPart. 

# Requirements #

In order to run SpecPart, the following requirements are mandatory: 
* [Julia](https://julialang.org/) version 1.6.0 or newer. 
* The following Julia libraries:
  * [Shuffle](https://docs.juliahub.com/Shuffle/X0eqg/0.1.1/)
  * [LightGraphs](https://github.com/sbromberger/LightGraphs.jl)
  * [SimpleWeightedGraphs](https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl)
  * [LDLFactorizations](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl)
  * [IterativeSolvers](https://iterativesolvers.julialinearalgebra.org/stable/)
  * [Laplacians](https://github.com/danspielman/Laplacians.jl)
  * [LinearMaps](https://github.com/JuliaLinearAlgebra/LinearMaps.jl)
* [CPLEX ILP solver](https://www.ibm.com/support/pages/downloading-ibm-ilog-cplex-optimization-studio-v1290) 
* [hMETIS](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview)
 
# SpecPart Parameters #

SpecPart accepts the following parameters:
* hg: The hypergraph file.
* pfile: The partition file.
* Nparts: The number of partitions. (Current implementation only allows bipartitions) 
* cycles: The number of random cycles used for graph construction. (Default 2)
* hyperedges_threshold: The threshold number of hyperedges used for running ILP. (Default 300)
* ub: The imbalance factor.
* nev: The number of eigenvectors to be computed. (Default 2)
* refine_iters: The number of ISSHP iterations. (Default 2)
* best_solns: The number of partition solutions picked for overlay-clustering. (Default 5)
* seed: The random seed for SpecPart.
* solver_iters: The number of iterations of LOBPCG. (Default 80)

# SpecPart Example #

```
julia> include("SpectralRefinement.jl")
julia> using Main.SpecPart
julia> SpecPart.SpectralRefinement(hg = "ibm01.hgr", pfile = "ibm01.hgr.2", Nparts = 2, cycles = 2, hyperedges_threshold = 300, ub = 2, nev = 2, refine_iters = 2, best_solns = 5)

[ Info: ================================================================================
[ Info: STARTING SUPERVISED SPECTRAL PARTITIONING ENGINE
[ Info: ================================================================================
[ Info: TOTAL VERTICES: 12752
[ Info: TOTAL HYPEREDGES: 14111
[ Info: POST PROCESSING :: TOTAL VERTICES: 12752
[ Info: POST PROCESSING :: TOTAL HYPEREDGES: 14111
[ Info: SIZE OF LARGEST HYPEREDGE: 42
[ Info: MAX CAPACITY CONSTRAINT: 6632
[ Info: MIN CAPACITY CONSTRAINT: 6120
[ Info: ================================================================================
[ Info: [HINT] CUT RECORDED IS 213
┌ Info: [EIGEN VECTOR DETAILS] :: Results of LOBPCG Algorithm
│  * Algorithm: LOBPCG - CholQR
│  * λ: [5.240705832350001e-8,2.0391573229227554e-7]
│  * Residual norm(s): [4.816542606411436e-16,6.027634032547232e-10]
│  * Convergence
│    * Iterations: 81
│    * Converged: false
└    * Iterations limit: 80
[ Info: [EIGEN VECTOR DETAILS] :: SOLVER TIME :: 3.860059666 seconds
[ Info: [SPECTRAL ITERATION 1] MIN CUT FOUND :: 330 :: TREE GEN-SOLVE TIME :: 10.548340258 seconds
┌ Info: [EIGEN VECTOR DETAILS] :: Results of LOBPCG Algorithm
│  * Algorithm: LOBPCG - CholQR
│  * λ: [5.074271058145485e-8,1.9977580740832276e-7]
│  * Residual norm(s): [6.266095133681997e-18,6.6504495135645066e-15]
│  * Convergence
│    * Iterations: 81
│    * Converged: false
└    * Iterations limit: 80
[ Info: [EIGEN VECTOR DETAILS] :: SOLVER TIME :: 0.860196429 seconds
[ Info: [SPECTRAL ITERATION 2] MIN CUT FOUND :: 328 :: TREE GEN-SOLVE TIME :: 0.474149347 seconds
[ Info: ================================================================================
[ Info: SIZE OF CLUSTERED HYPERGRAPH IS: 430 VERTICES AND 636 HYPEREDGES
[ Info: ================================================================================
[ Info: PARTITIONING CLUSTERED HYPERGRAPH .....
[ Info: RUNNING CPLEX AS GOLDEN PARTITIONER
[ Info: [POST TOOL CUT] CUT RECORDED IS 254
[ Info: ================================================================================
[ Info: [SUPERVISED SPECTRAL] CUTSIZE OBTAINED: 213
[ Info: [SUPERVISED SPECTRAL] AREA SPLIT OBTAINED: 6500 [0.51%] :: 6252 [0.49%]
[ Info: ================================================================================
[ Info: [RUNTIME] IO PROCESSING :: 0.297085323 SECONDS
[ Info: [RUNTIME] SPECTRAL :: 16.228560581 SECONDS
[ Info: [RUNTIME] CLUSTERING :: 0.157013099 SECONDS
[ Info: [RUNTIME] PARTITION CLUSTERED HYPERGRAPH :: 1.082736077 SECONDS
[ Info: ================================================================================
[ Info: [RUNTIME] TOTAL EXECUTION TIME :: 17.76539508 SECONDS
[ Info: ================================================================================
```

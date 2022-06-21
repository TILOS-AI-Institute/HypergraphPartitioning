# SpecPart User Guide #

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

Libraries can be added by doing the following in Julia REPL:

```
julia> using Pkg
julia> Pkg.add("Package name")
```

# SpecPart Parameters #

SpecPart accepts the following parameters:

| Parameter   | Description |
| ----------- | ----------- |
| hg      | The hypergraph file.       |
| pfile   | The partition file.        |
| Nparts   | The number of partitions. (Current implementation only allows bipartitions)         |
| cycles   | The number of random cycles used for graph construction. (Default 2)        |
| hyperedges_threshold   | The threshold number of hyperedges used for running ILP. (Default 300)        |
| ub   | The imbalance factor.        |
| nev   | The number of eigenvectors to be computed. (Default 2)        |
| refine_iters   | The number of ISSHP iterations. (Default 2)        |
| best_solns   | The number of partition solutions picked for overlay-clustering. (Default 5)        |
| seed   | The random seed for SpecPart.        |
| solver_iters   | The number of iterations of LOBPCG. (Default 80)        |

# SpecPart Example #
We show how we run SpecPart on ISPD98 testcase "ibm02.hgr" with an input partition. As seen from this example SpecPart improves the initial partition of 339 to 334.

```
julia> include("SpecPart/SpectralRefinement.jl")
julia> using Main.SpecPart
julia> SpecPart.SpectralRefinement(hg = "benchmark/ISPD_benchmark/ibm02.hgr", pfile = "solutions/ISPD_benchmark_solutions/hMetis/UBfactor_2/ibm02.hgr.k.2.UBfactor.2.seed.0", 
Nparts = 2, cycles = 2, hyperedges_threshold = 300, ub = 2, nev = 2, refine_iters = 2, best_solns = 5)

[ Info: ================================================================================
[ Info: STARTING SPECPART: SUPERVISED SPECTRAL PARTITIONING ENGINE
[ Info: ================================================================================
[ Info: TOTAL VERTICES: 19601
[ Info: TOTAL HYPEREDGES: 19584
[ Info: POST PROCESSING :: TOTAL VERTICES: 19601
[ Info: POST PROCESSING :: TOTAL HYPEREDGES: 19584
[ Info: SIZE OF LARGEST HYPEREDGE: 134
[ Info: MAX CAPACITY CONSTRAINT: 10193
[ Info: MIN CAPACITY CONSTRAINT: 9408
[ Info: ================================================================================
[ Info: [HINT] CUT RECORDED IS 339
┌ Info: [EIGEN VECTOR DETAILS] :: Results of LOBPCG Algorithm
│  * Algorithm: LOBPCG - CholQR
│  * λ: [6.09867746662552e-8,8.185494183481112e-8]
│  * Residual norm(s): [1.3833209684774178e-17,1.3317377476442714e-16]
│  * Convergence
│    * Iterations: 81
│    * Converged: false
└    * Iterations limit: 80
[ Info: [EIGEN VECTOR DETAILS] :: SOLVER TIME :: 1.712944723 seconds
[ Info: [SPECTRAL ITERATION 1] MIN CUT FOUND :: 434 :: TREE GEN-SOLVE TIME :: 1.092160598 seconds
┌ Info: [EIGEN VECTOR DETAILS] :: Results of LOBPCG Algorithm
│  * Algorithm: LOBPCG - CholQR
│  * λ: [6.063794065520354e-8,8.233442427641748e-8]
│  * Residual norm(s): [9.10638903776076e-16,1.0476424071977864e-14]
│  * Convergence
│    * Iterations: 81
│    * Converged: false
└    * Iterations limit: 80
[ Info: [EIGEN VECTOR DETAILS] :: SOLVER TIME :: 1.905182079 seconds
[ Info: [SPECTRAL ITERATION 2] MIN CUT FOUND :: 431 :: TREE GEN-SOLVE TIME :: 0.68814692 seconds
[ Info: ================================================================================
[ Info: SIZE OF CLUSTERED HYPERGRAPH IS: 952 VERTICES AND 875 HYPEREDGES
[ Info: ================================================================================
[ Info: PARTITIONING CLUSTERED HYPERGRAPH .....
[ Info: RUNNING CPLEX AS GOLDEN PARTITIONER
[ Info: [POST TOOL CUT] CUT RECORDED IS 334
[ Info: ================================================================================
[ Info: [SUPERVISED SPECTRAL] CUTSIZE OBTAINED: 334
[ Info: [SUPERVISED SPECTRAL] AREA SPLIT OBTAINED: 10191 [0.52%] :: 9410 [0.48%]
[ Info: ================================================================================
[ Info: [RUNTIME] IO PROCESSING :: 0.11667171 SECONDS
[ Info: [RUNTIME] SPECTRAL :: 5.562175724 SECONDS
[ Info: [RUNTIME] CLUSTERING :: 0.048536965 SECONDS
[ Info: [RUNTIME] PARTITION CLUSTERED HYPERGRAPH :: 1.311114935 SECONDS
[ Info: ================================================================================
[ Info: [RUNTIME] TOTAL EXECUTION TIME :: 7.038499334000001 SECONDS
[ Info: ================================================================================

```

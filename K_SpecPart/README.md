# K-SpecPart #

K-SPecPart is an extension (and improvement) over SpecPart that can compute multi-way (K >= 2) partitions. 

# Requirements #

In order to run K-SpecPart, the following requirements are mandatory: 
* [Julia](https://julialang.org/) version 1.6.0 or newer. 
* The following Julia libraries:
  * [Shuffle](https://docs.juliahub.com/Shuffle/X0eqg/0.1.1/)
  * [LightGraphs](https://github.com/sbromberger/LightGraphs.jl)
  * [SimpleWeightedGraphs](https://github.com/JuliaGraphs/SimpleWeightedGraphs.jl)
  * [LDLFactorizations](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl)
  * [IterativeSolvers](https://iterativesolvers.julialinearalgebra.org/stable/)
  * [Laplacians](https://github.com/danspielman/Laplacians.jl)
  * [LinearMaps](https://github.com/JuliaLinearAlgebra/LinearMaps.jl)
* [OR-Tools](https://developers.google.com/optimization)
* If you have acces to [CPLEX ILP solver](https://www.ibm.com/support/pages/downloading-ibm-ilog-cplex-optimization-studio-v1290), you can use CPLEX in place of OR-Tools. 
* [METIS](https://github.com/KarypisLab/METIS)
* [hMETIS](http://glaros.dtc.umn.edu/gkhome/metis/hmetis/overview)
* [FM](https://github.com/ABKGroup/TritonPart_OpenROAD.git)

Julia libraries can be installed by doing the following in Julia REPL:

```
julia> using Pkg
julia> Pkg.add("Package name")
```

# K-SpecPart Parameters #

K-SpecPart accepts the following parameters:

| Parameter   | Description |
| ----------- | ----------- |
| hypergraph_fixed_file      | The fixed file consisting of fixed vertices       |
| hint_file   | The hint partition file        |
| imb  | The imbalance factor         |
| num_parts   | The number of partitions        |
| eigvecs   | The number of eigenvectors (Default 2)        |
| refine_iters   | The number of K-SpecPart iterations (Default 2)        |
| solver_iters   | The number of LOCPCG solver iterations (Default 40)        |
| best_solns   | The number of partition solutions picked for overlay-clustering (Default 5)        |
| ncycles   | The number of random cycles for constructing spectral sparsifier graph (Default 2)        |
| seed   | The random seed for K-SpecPart (Default 0)       |
